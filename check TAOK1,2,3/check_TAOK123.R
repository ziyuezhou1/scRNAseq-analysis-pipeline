# Annotation
# Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(scran)
library(Seurat)
library(SingleCellExperiment)

# Loading data
so <- readRDS(file.path("./", "so.renamed_res0.2.rds"))
pdf("./Taok1_split_group_vln.pdf",
    width = 40,
    height = 8)
VlnPlot(object = so, features = c("Taok1"), split.by = "group_id", ncol = 1, pt.size = 0)
dev.off()
pdf("./Taok1_split_group_dim.pdf",
    width = 40,
    height = 8)
FeaturePlot(so, features = c("Taok1"), reduction = "umap", split.by = "group_id", coord.fixed = 1)
dev.off()
pdf("./Taok2_split_group_vln.pdf",
    width = 40,
    height = 8)
VlnPlot(object = so, features = c("Taok2"), split.by = "group_id", ncol = 1, pt.size = 0)
dev.off()
pdf("./Taok3_split_group_vln.pdf",
    width = 40,
    height = 8)
VlnPlot(object = so, features = c("Taok3"), split.by = "group_id", ncol = 1, pt.size = 0)
dev.off()
