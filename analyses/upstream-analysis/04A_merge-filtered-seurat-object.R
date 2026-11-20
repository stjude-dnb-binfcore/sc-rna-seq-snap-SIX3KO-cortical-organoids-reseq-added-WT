# 09/04/2026

# This script is to merge non-GEMx and GEMx samples.
# The filtered object generated from `/sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq/sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq-gene-300/analyses/upstream-analysis/results/04_Filter_object/merge-filtered-seurat-obj.R`
# was used to merge with the filtered object from `/sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq-added-WT/analyses/upstream-analysis/results/04_Filter_object/seurat_obj_merged_filtered.rds`

setwd("/sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq-added-WT/")
library(Seurat)


added_wt <- readRDS("./analyses/upstream-analysis/results/04_Filter_object/seurat_obj_merged_filtered-gene-1000-wo-outlier-threshold.rds")
prev_obj <- readRDS("/sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq/sc-rna-seq-snap-SIX3KO-cortical-organoids-reseq-gene-300/analyses/upstream-analysis/results/04_Filter_object/seurat_obj_merged_filtered.rds")

seurat_obj_merged_filtered_cust <- merge(x = added_wt, y = prev_obj)


DefaultAssay(seurat_obj_merged_filtered_cust) <- "RNA"

seurat_obj_merged_filtered_cust <- NormalizeData(seurat_obj_merged_filtered_cust)
seurat_obj_merged_filtered_cust <- FindVariableFeatures(seurat_obj_merged_filtered_cust)

seurat_obj_merged_filtered_cust <- ScaleData(seurat_obj_merged_filtered_cust, features = rownames(seurat_obj_merged_filtered_cust))

seurat_obj_merged_filtered_cust <- RunPCA(seurat_obj_merged_filtered_cust,
                                          features = VariableFeatures(seurat_obj_merged_filtered_cust))

saveRDS(seurat_obj_merged_filtered_cust, "./analyses/upstream-analysis/results/04_Filter_object/seurat_obj_merged_filtered.rds")
