# loading packages
library(Seurat)

# loading data
seu <- readRDS("../../data/so.renamed_res0.2.rds")

VlnPlot(seu, features = c("Tnf"), split.by = "group_id")
VlnPlot(seu, features = c("Il6"), split.by = "group_id")
VlnPlot(seu, features = c("Il1a"), split.by = "group_id")
VlnPlot(seu, features = c("Il1b"), split.by = "group_id")
VlnPlot(seu, features = c("Ccl2"), split.by = "group_id")
VlnPlot(seu, features = c("Ccl5"), split.by = "group_id")
VlnPlot(seu, features = c("Cxcl1"), split.by = "group_id")
VlnPlot(seu, features = c("Il12a"), split.by = "group_id")
