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
# write.table(exprMat, file = "./exprMat.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
cellInfo <- seu@meta.data
add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }

  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")

  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }

  invisible(loom)
}

loom <- build_loom("so.renamed_res.0.2.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

# download ref files
# if (!dir.exists("cisTarget_databases")) {
#   dir.create("cisTarget_databases"); setwd("cisTarget_databases")
# }
# options(timeout = 6000)
# # dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
# #              "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
# #              "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
# #              "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather")
# dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
#              "https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# # mc9nr: Motif collection version 9: 24k motifs
# for(featherURL in dbFiles)
# {
#   download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
# }
### Load data
# loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")

loom <- open_loom("so.renamed_res.0.2.loom")
exprMat <- get_dgem(loom)
# write.csv(t(exprMat), file = "./exprMat.csv", row.names = TRUE, col.names = TRUE)
write.table(exprMat, file = "./exprMat.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
# cellInfo <- get_cell_annotation(loom)
# close_loom(loom)
# 
# library(sceasy)
# library(reticulate)
# loompy <- reticulate::import('loompy')
# seu <- readRDS(file = "so.renamed_res0.2.rds")
# sceasy::convertFormat(obj = seu, from = "seurat", to = "sce", outFile = "sce.renamed_res.0.2.rds")
# sceasy::convertFormat(obj = sce, from = "sce", to = "loom", outFile = "so.renamed_res.0.2.loom")
