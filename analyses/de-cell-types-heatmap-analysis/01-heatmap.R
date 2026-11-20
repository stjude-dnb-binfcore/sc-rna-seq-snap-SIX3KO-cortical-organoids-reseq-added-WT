
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)


# Read the gene-list data frame
gene_lists <- read.delim("genelist-of-interest.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
seurat_obj <- readRDS("analyses/cell-types-annotation/results/cell_types_annotations_all/seurat_obj_cell_types_annotations_all.rds")

# Rename Glial progenitor directly in the Seurat metadata
seurat_obj$curated.cell.types <- as.character(seurat_obj$curated.cell.types)
seurat_obj$curated.cell.types[
  seurat_obj$curated.cell.types == "Glial progenitor"] <- "Mesenchyme Cells"

# Rename the corresponding gene-list column if needed
if ("Glial progenitor" %in% colnames(gene_lists)) {
  colnames(gene_lists)[colnames(gene_lists) == "Glial progenitor"] <- "Mesenchyme Cells"
}

# Cell types to include and desired order
celltypes_keep <- c("Neuroepithelium","Subcortical neuron",
                    "Intermediate progenitor", "Mesenchyme Cells")

# Convert gene lists to long format
gene_map <- gene_lists %>%
  select(all_of(celltypes_keep)) %>%
  pivot_longer(cols = everything(), names_to = "marker_cell_type", values_to = "gene") %>%
  filter(!is.na(gene), gene != "") %>% distinct(marker_cell_type, gene)

# Keep genes present in Seurat object
genes <- unique(gene_map$gene)
genes <- intersect(genes, rownames(seurat_obj))

gene_map <- gene_map %>% filter(gene %in% genes)

# Extract normalized expression + metadata
expr_df <- FetchData(seurat_obj,
  vars = c(genes, "orig.ident", "condition", "curated.cell.types"), slot = "data")

# Keep selected cell types and convert to long format
expr_long <- expr_df %>% rownames_to_column("cell") %>%
  filter(curated.cell.types %in% celltypes_keep) %>%
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression")

# Average normalized expression per sample x cell type
expr_summary <- expr_long %>%
  group_by(curated.cell.types, orig.ident, condition, gene) %>%
  summarise(avg_expression = mean(expression),
    .groups = "drop")

# Scale each gene across ALL sample x cell-type combinations
expr_summary <- expr_summary %>%
  group_by(gene) %>% mutate(scaled_expression = as.numeric(scale(avg_expression))) %>%
  ungroup()

# Create unique heatmap column ID
expr_summary <- expr_summary %>%
  mutate(group_id = paste(curated.cell.types, orig.ident, sep = " | "))

# Create heatmap matrix
heatmap_mat <- expr_summary %>%
  select(gene, group_id, scaled_expression) %>%
  pivot_wider(names_from = group_id, values_from = scaled_expression) %>%
  column_to_rownames("gene") %>% as.matrix()

# Create column annotation information
annotation_col <- expr_summary %>% select(group_id, orig.ident, condition, curated.cell.types) %>%
  distinct() %>% as.data.frame()

# Add short sample names
annotation_col <- annotation_col %>%
  mutate(sample_alt = case_when(
      orig.ident == "Six3_D28_CtrlC" ~ "C",
      orig.ident == "Six3_D28_CtrlD" ~ "D",
      orig.ident == "Six3_CtrlD_1" ~ "D1",
      orig.ident == "Six3_CtrlD_2" ~ "D2",
      orig.ident == "Six3_CtrlD_3" ~ "D3",
      orig.ident == "Six3_D28_125_P" ~ "125P",
      orig.ident == "Six3_D28_125_X" ~ "125X",
      orig.ident == "Six3_D28_24_P" ~ "24P",
      orig.ident == "Six3_D28_76M" ~ "76M",
      orig.ident == "Six3_D28_209" ~ "209",
      TRUE ~ orig.ident))

# Desired condition and cell-type order
annotation_col$condition <- factor(annotation_col$condition, levels = c("WT", "KO"))

annotation_col$curated.cell.types <- factor(annotation_col$curated.cell.types, levels = celltypes_keep)

# Desired sample order
sample_order_levels <- c("C", "D", "D1", "D2", "D3",
  "125P", "125X", "24P", "76M", "209")

annotation_col$sample_alt <- factor(annotation_col$sample_alt, levels = sample_order_levels)

# Order columns:
# Cell type -> WT/KO -> desired sample order
annotation_col <- annotation_col %>% arrange(curated.cell.types, condition, sample_alt)

column_order <- annotation_col$group_id

# Reorder heatmap columns
heatmap_mat <- heatmap_mat[, column_order, drop = FALSE]

# Use unique group IDs as row names for annotation table
rownames(annotation_col) <- annotation_col$group_id
annotation_col <- annotation_col[column_order, , drop = FALSE]

# Visible heatmap column labels = short sample names
labels_col <- as.character(annotation_col$sample_alt)

# Keep only annotation bars wanted in plot
annotation_col <- annotation_col %>% select(condition, curated.cell.types)

# Rename column annotation field
colnames(annotation_col)[colnames(annotation_col) == "curated.cell.types"] <- "Cell type"

# Order genes by marker cell type
gene_map$marker_cell_type <- factor(gene_map$marker_cell_type, levels = celltypes_keep)

gene_order <- gene_map %>% arrange(marker_cell_type) %>% distinct(gene, .keep_all = TRUE) %>% pull(gene)

gene_order <- gene_order[gene_order %in% rownames(heatmap_mat)]

# Reorder heatmap rows
heatmap_mat <- heatmap_mat[gene_order, , drop = FALSE]

# Create row annotation
annotation_row <- gene_map %>%
  arrange(marker_cell_type) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(gene %in% rownames(heatmap_mat)) %>%
  select(gene, marker_cell_type) %>%
  as.data.frame()

rownames(annotation_row) <- annotation_row$gene
annotation_row$gene <- NULL

# Match heatmap row order
annotation_row <- annotation_row[rownames(heatmap_mat), , drop = FALSE]

# Rename row annotation field
colnames(annotation_row)[colnames(annotation_row) == "marker_cell_type"] <- "Cell type"

# Add gaps between cell-type column blocks
celltype_groups_col <- as.character(annotation_col$`Cell type`)

gaps_col <- cumsum(rle(celltype_groups_col)$lengths)

gaps_col <- gaps_col[-length(gaps_col)]

# Add gaps between marker-gene row blocks
celltype_groups_row <- as.character(annotation_row$`Cell type`)

gaps_row <- cumsum(
  rle(celltype_groups_row)$lengths)

gaps_row <- gaps_row[
  -length(gaps_row)]

# Cell-type colors used for BOTH row and column annotations
celltype_colors <- c(
  "Neuroepithelium" = "#B0D2DF",
  "Subcortical neuron" = "#E8BE3E",
  "Intermediate progenitor" = "#6A329F",
  "Mesenchyme Cells" = "#B4A7D6")

# Annotation colors
annotation_colors <- list(
  condition = c(
    "WT" = "#2C7BB6",
    "KO" = "#D7191C"),
  `Cell type` = celltype_colors)


# Save as SVG
svg(filename = "analyses/custom-figures/heatmap/Scaled_average_normalized_expression.svg",
  width = 28,
  height = 14)

pheatmap(heatmap_mat,
  scale = "none",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  labels_col = labels_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  border_color = NA,
  cellwidth = 18,
  cellheight = 8,
  fontsize_row = 7,
  fontsize_col = 9,
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  angle_col = 45,
  main = "Gene expression across samples and cell types")

dev.off()
