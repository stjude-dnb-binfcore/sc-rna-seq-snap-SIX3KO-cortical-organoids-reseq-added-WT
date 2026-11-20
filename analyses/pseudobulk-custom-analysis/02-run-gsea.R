library(dplyr)
library(readr)
library(msigdbr)
library(fgsea)

setwd("sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq-added-WT/analyses/")

# Input/output directories
pseudobulk_dir <- "cluster-cell-calling/results/pseudo-bulk"
gsea_dir <- file.path(pseudobulk_dir, "GSEA")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

# Find pseudobulk DE TSV files
de_files <- list.files(
  path = pseudobulk_dir,
  pattern = "^Pseudobulk_DE_.*\\.tsv$",
  full.names = TRUE
)

if (length(de_files) == 0) {
  stop("No pseudobulk DE TSV files found in: ", pseudobulk_dir)
}

cat("Found", length(de_files), "pseudobulk DE files\n")
print(de_files)

# Gene sets
Species <- "Homo sapiens"

hallmark_df <- msigdbr(species = Species, category = "H")
reactome_df <- msigdbr(species = Species, category = "C2", subcategory = "CP:REACTOME")
go_bp_df <- msigdbr(species = Species, category = "C5", subcategory = "GO:BP")

hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
reactome_sets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
go_bp_sets <- split(go_bp_df$gene_symbol, go_bp_df$gs_name)

cat("Hallmark gene sets:", length(hallmark_sets), "\n")
cat("Reactome gene sets:", length(reactome_sets), "\n")
cat("GO Biological Process gene sets:", length(go_bp_sets), "\n")

# Store results
hallmark_results <- list()
reactome_results <- list()
go_bp_results <- list()

# Run GSEA for each cell-type DE file
for (de_file in de_files) {
  
  de_name <- tools::file_path_sans_ext(basename(de_file))
  de_name <- gsub(" ", "_", de_name)
  
  cat("\nRunning GSEA:", de_name, "\n")
  
  de <- read_tsv(de_file, show_col_types = FALSE)
  
  if (!all(c("Gene", "avg_log2FC") %in% colnames(de))) {
    warning("Skipping ", de_name, ": Gene or avg_log2FC column missing")
    next
  }
  
  # Use all genes with valid avg_log2FC
  de_ranked <- de %>%
    filter(!is.na(Gene), !is.na(avg_log2FC))
  
  ranks <- de_ranked$avg_log2FC
  names(ranks) <- de_ranked$Gene
  
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  
  cat("Genes used for GSEA:", length(ranks), "\n")
  
  # Hallmark
  hallmark_gsea <- fgsea(pathways = hallmark_sets,
    stats = ranks,
    minSize = 15,
    maxSize = 500) %>%
    as.data.frame() %>%
    arrange(padj)
  
  hallmark_gsea$leadingEdge <- vapply(hallmark_gsea$leadingEdge,
    function(x) paste(x, collapse = ","),
    character(1))
  
  hallmark_results[[de_name]] <- hallmark_gsea
  
  write.table(hallmark_gsea,
    file = file.path(gsea_dir, paste0(de_name, "_Hallmark_GSEA.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)
  
  # Reactome
  reactome_gsea <- fgsea(
    pathways = reactome_sets,
    stats = ranks,
    minSize = 15,
    maxSize = 500) %>%
    as.data.frame() %>%
    arrange(padj)
  
  reactome_gsea$leadingEdge <- vapply(
    reactome_gsea$leadingEdge,
    function(x) paste(x, collapse = ","),
    character(1))
  
  reactome_results[[de_name]] <- reactome_gsea
  
  write.table(
    reactome_gsea,
    file = file.path(gsea_dir, paste0(de_name, "_Reactome_GSEA.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)
  
  # GO Biological Process
  go_bp_gsea <- fgsea(
    pathways = go_bp_sets,
    stats = ranks,
    minSize = 15,
    maxSize = 500) %>%
    as.data.frame() %>%
    arrange(padj)
  
  go_bp_gsea$leadingEdge <- vapply(
    go_bp_gsea$leadingEdge,
    function(x) paste(x, collapse = ","),
    character(1))
  
  go_bp_results[[de_name]] <- go_bp_gsea
  
  write.table(
    go_bp_gsea,
    file = file.path(gsea_dir, paste0(de_name, "_GO_BP_GSEA.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)
  
  # Print results
  cat("\nHallmark:\n")
  print(hallmark_gsea %>% select(pathway, NES, pval, padj, size, leadingEdge))
  
  cat("\nReactome:\n")
  print(reactome_gsea %>% select(pathway, NES, pval, padj, size, leadingEdge))
  
  cat("\nGO Biological Process:\n")
  print(go_bp_gsea %>% select(pathway, NES, pval, padj, size, leadingEdge))
}

# Combined Hallmark
all_hallmark <- bind_rows(hallmark_results, .id = "DE_comparison")
write.table(
  all_hallmark,
  file = file.path(gsea_dir, "ALL_celltypes_Hallmark_GSEA.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

# Combined Reactome
all_reactome <- bind_rows(reactome_results, .id = "DE_comparison")
write.table(
  all_reactome,
  file = file.path(gsea_dir, "ALL_celltypes_Reactome_GSEA.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

# Combined GO Biological Process
all_go_bp <- bind_rows(go_bp_results, .id = "DE_comparison")
write.table(
  all_go_bp,
  file = file.path(gsea_dir, "ALL_celltypes_GO_BP_GSEA.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

# Significant combined results
all_hallmark_sig <- all_hallmark %>% filter(!is.na(padj), padj < 0.05)
all_reactome_sig <- all_reactome %>% filter(!is.na(padj), padj < 0.05)
all_go_bp_sig <- all_go_bp %>% filter(!is.na(padj), padj < 0.05)

write.table(
  all_hallmark_sig,
  file = file.path(gsea_dir, "ALL_celltypes_Hallmark_GSEA_significant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

write.table(
  all_reactome_sig,
  file = file.path(gsea_dir, "ALL_celltypes_Reactome_GSEA_significant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

write.table(
  all_go_bp_sig,
  file = file.path(gsea_dir, "ALL_celltypes_GO_BP_GSEA_significant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

# Print summary
cat("\nGSEA complete\n")
cat("Cell-type DE files processed:", length(hallmark_results), "\n")
cat("Results saved to:", gsea_dir, "\n")
