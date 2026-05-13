library(monocle)
library(Seurat)
# load seurat object
seu <- readRDS("../../#5-DEG - no FC filtered/data/so.renamed_res0.2.rds")
groups <- c("ZT6", "ZT18","SD3D 8h/day (ZT0-ZT8)")
seu <- subset(x = seu, idents = c("Inh-VIP&Lamp5","Interneuron"), subset = group_id %in% groups, invert = FALSE)

# Idents(seu) <- 'cluster'
DefaultAssay(seu) <- 'RNA'

# transfer into cds
cds <- as.CellDataSet(seu)

## Monocle2 process
# estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
