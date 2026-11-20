# Pseudobulk DE analysis on curated clusters - WT vs SIX3
library(Seurat)
library(DESeq2)
library(tibble)
library(Matrix)
library(dplyr)

setwd("sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq-added-WT/analyses/cluster-cell-calling/results/pseudo-bulk/")

seurat_obj = readRDS("sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq-added-WT/analyses/cell-types-annotation/results/cell_types_annotations_all/seurat_obj_cell_types_annotations_all.rds")

unique(seurat_obj$ID)
unique(seurat_obj$curated.cell.types)
unique(seurat_obj$condition)

Generate.Pseudobulk.DE <- function(Seurat_Obj, Sample.Column = "ID", Subset.Column = "curated.cell.types",
                                   Idents.Column = "condition", Ident.1 = "KO", Ident.2 = "WT", Results.Path = getwd()) 
  {
  # Create pseudobulk ID
  Seurat_Obj$PseudoBulk_ID <- paste(Seurat_Obj@meta.data[[Sample.Column]], Seurat_Obj@meta.data[[Subset.Column]], Seurat_Obj@meta.data[[Idents.Column]],
    sep = "_")
  
  # Keep metadata for each pseudobulk sample
  PB_Metadata <- unique(Seurat_Obj@meta.data[, c("PseudoBulk_ID", Sample.Column,
                                                 Subset.Column, Idents.Column)])
  
  rownames(PB_Metadata) <- PB_Metadata$PseudoBulk_ID
  
  # Aggregate raw counts
  Pseudo_Obj <- AggregateExpression(Seurat_Obj, return.seurat = TRUE,
                                    assays = "RNA",slot =  "counts", group.by = "PseudoBulk_ID")
  
  
  # Add metadata back to pseudobulk object
  PB_IDs <- Cells(Pseudo_Obj)
  
  Pseudo_Obj[[Sample.Column]] <- PB_Metadata[PB_IDs, Sample.Column]
  
  Pseudo_Obj[[Subset.Column]] <- PB_Metadata[PB_IDs, Subset.Column]
  
  Pseudo_Obj[[Idents.Column]] <- PB_Metadata[PB_IDs, Idents.Column]
  
  
  # Column to iterate over
  Pseudo_Obj[["Subset.Variable"]] <- Pseudo_Obj[[Subset.Column]]
  
  # Run DE separately for each curated cell type
  for (i in unique(Pseudo_Obj$Subset.Variable)) {
    
    print(i)
    
    Bulk_Cluster <- subset( Pseudo_Obj, Subset.Variable == i)
    
    Idents(Bulk_Cluster) <- Idents.Column
    
    # Number of biological pseudobulk samples
    print(table(Idents(Bulk_Cluster)))
    
    # KO vs WT
    DE_Genes <- FindMarkers(Bulk_Cluster, ident.1 = Ident.1, ident.2 = Ident.2,
                            slot = "counts", test.use = "DESeq2", min.cells.group = 2)
    
    # Add gene names
    DE_Genes$Gene <- rownames(DE_Genes)
    
    DE_Genes <- DE_Genes[, c("Gene", "p_val", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")]
    
    
    # Write DE results
    write.table(DE_Genes, file = file.path(Results.Path, paste0("Pseudobulk_DE_", i, "_", Ident.1, "_vs_", Ident.2, ".tsv")), 
                quote = FALSE, row.names = FALSE,col.names = TRUE, sep = "\t")
  }
  
  return(Pseudo_Obj)
}



Pseudo_Obj <- Generate.Pseudobulk.DE(Seurat_Obj = seurat_obj, Results.Path = getwd())


######################################################################################

# Additional aanlysis: Per-sample pseudobulk log2 fold-change analysis

sample_col <- "orig.ident"
condition_col <- "condition"
celltype_col <- "curated.cell.types"
top_n <- 20

output_dir <- "pseudobulk_replicate_expression"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get raw counts and metadata
counts <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
meta <- seurat_obj@meta.data

# Create sample x cell-type labels
meta$sample_celltype <- paste(meta[[sample_col]], meta[[celltype_col]], sep = "__")

# Aggregate raw counts by sample x cell type
groups <- unique(meta$sample_celltype)

pseudobulk_counts <- sapply(groups, function(grp) {
  cells <- rownames(meta)[meta$sample_celltype == grp]
  Matrix::rowSums(counts[, cells, drop = FALSE])
})

pseudobulk_counts <- as.matrix(pseudobulk_counts)
colnames(pseudobulk_counts) <- groups

# Create pseudobulk metadata
pb_meta <- meta %>%
  dplyr::select(all_of(sample_col), all_of(condition_col),
                all_of(celltype_col), sample_celltype) %>%
  dplyr::distinct() %>%
  as.data.frame()

rownames(pb_meta) <- pb_meta$sample_celltype
pb_meta <- pb_meta[colnames(pseudobulk_counts), , drop = FALSE]

# Cell numbers per sample x cell type
cell_number_table <- meta %>%
  dplyr::count(.data[[condition_col]], .data[[sample_col]],
               .data[[celltype_col]], name = "n_cells")

write.table(
  cell_number_table,
  file.path(output_dir, "cell_numbers_per_sample_celltype.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Run for all cell types
celltypes <- sort(unique(pb_meta[[celltype_col]]))

all_DE_results <- list()
all_normalized_expression <- list()

for (ct in celltypes) {
  
  cat("Processing:", ct, "\n")
  
  keep <- pb_meta[[celltype_col]] == ct
  counts_ct <- pseudobulk_counts[, keep, drop = FALSE]
  meta_ct <- pb_meta[keep, , drop = FALSE]
  
  condition_counts <- table(meta_ct[[condition_col]])
  print(condition_counts)
  
  if (!all(c("KO", "WT") %in% names(condition_counts))) {
    warning("Skipping ", ct, ": KO and WT are not both represented.")
    next
  }
  
  # Keep genes with >=10 raw counts in at least 2 pseudobulk samples
  keep_genes <- rowSums(counts_ct >= 10) >= 2
  counts_ct <- counts_ct[keep_genes, , drop = FALSE]
  
  # DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_ct),
    colData = meta_ct,
    design = ~ condition
  )
  
  dds$condition <- relevel(factor(dds$condition), ref = "WT")
  dds <- DESeq2::DESeq(dds)
  
  res <- DESeq2::results(dds, contrast = c("condition", "KO", "WT"))
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::arrange(padj)
  
  res_df$CellType <- ct
  res_df <- res_df %>% dplyr::select(CellType, Gene, everything())
  
  all_DE_results[[ct]] <- res_df
  
  ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  write.table(
    res_df,
    file.path(output_dir, paste0(ct_safe, "_KO_vs_WT_DESeq2.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # DESeq2 normalized pseudobulk expression for ALL genes
  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  colnames(normalized_counts) <- paste(
    meta_ct[[sample_col]],
    gsub("[^A-Za-z0-9]+", "_", ct),
    meta_ct[[condition_col]],
    sep = "_"
  )
  
  normalized_df <- data.frame(
    Gene = rownames(normalized_counts),
    normalized_counts,
    check.names = FALSE
  )
  
  normalized_df$CellType <- ct
  normalized_df <- normalized_df %>% dplyr::select(CellType, Gene, everything())
  rownames(normalized_df) <- NULL
  
  all_normalized_expression[[ct]] <- normalized_df
  
  write.table(
    normalized_df,
    file.path(output_dir, paste0(ct_safe, "_ALL_genes_normalized_expression.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

