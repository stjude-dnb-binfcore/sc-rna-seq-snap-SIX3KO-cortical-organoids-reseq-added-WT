# Generate Box Plot/Violin Plot

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

seurat_obj =  readRDS("/cell-types-annotation/results-gene-1000-wo-outlier-threshold-curatedlabels-added/cell_types_annotations_all/seurat_obj_cell_types_annotations_all.rds")
seurat_obj$curated.cell.types <- as.character(seurat_obj$curated.cell.types)

seurat_obj$curated.cell.types[seurat_obj$curated.cell.types == "Glial progenitor"] <- "Mesenchyme Cells"

# Get mesenchymal cells only
mesenchyme <- subset(seurat_obj, subset = curated.cell.types == "Mesenchyme Cells")

# Replicate-level Version
mesenchyme_meta <- mesenchyme_meta %>%
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

mesenchyme_meta$sample_alt <- factor(
  mesenchyme_meta$sample_alt,
  levels = c(
    "C",
    "D",
    "D1",
    "D2",
    "D3",
    "125X",
    "125P",
    "24P",
    "76M",
    "209"))


# Genes detected by sample
p4 <- ggplot(mesenchyme_meta,
  aes(x = sample_alt, y = nFeature_RNA, fill = condition)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  labs(
    x = "Sample",
    y = "Genes detected per nucleus",
    title = "Genes detected by sample") +
  theme_classic()

sample_plot <- p4

svg(filename = "Mesenchymal-SIX3KO-vs-WT/Mesenchymal_QC_by_sample.svg",
  width = 14,
  height = 6)

print(sample_plot)

dev.off()
