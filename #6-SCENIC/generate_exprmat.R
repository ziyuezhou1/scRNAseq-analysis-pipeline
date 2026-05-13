# loading packages
library(loomR)
library(Seurat)
library(dplyr)
# library(scater)
library(hdf5r)
library(iterators)
# library(ggplot2)
library(cowplot)
library(Matrix)
library(SeuratDisk)
library(SCopeLoomR)

# Convert rds to loom
seu <- readRDS(file = "so.renamed_res0.2.rds")
# table(Idents(seu))
# write.csv(t(as.matrix(seu@assays$RNA@counts)), file = "exprMat.csv")
exprMat <- seu@assays$RNA@data
write.table(exprMat, file = "./exprMat.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)