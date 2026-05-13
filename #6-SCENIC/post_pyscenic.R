# do post-pyscenic analysis

# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(Seurat)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

### Define study-specific parameters

out_dir <- "./outs_2023-04-12"
out_files <- list.files(out_dir, full.names = TRUE)

seurat_path <- "so.renamed_res0.2.rds"
analysis <- "so.renamed_res0.2_non_integrated"


### Load pySCENIC output loom file (regulons)

# from pySCENIC working directory:

path2loom <- out_files[grep("pySCENIC_viz.loom", out_files)]
# using *filtered* loom to extract regulons
loom <- open_loom(path2loom, mode = "r")

# Read pySCENIC information from loom file
# WARNING: column attribute names may be different depending on the SCopeLoomR/pySCENIC version used!
exprMat <- get_dgem(loom)
# exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
# => regulons
regulons_incidMat <- get_regulons(loom,  # as incidence matrix
                            column.attr.name = "Regulons")
dim(regulons_incidMat)
write.csv(regulons_incidMat,
          file = paste0(out_dir, "/", analysis, "_pySCENIC_regulons_incidence_matrix.csv"))

regulons <- regulonsToGeneLists(regulons_incidMat)  # can convert to gene list for processing with R
head(regulons)

# Regulon AUC and thresholds

regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)

# Embeddings (tsne/umap)
embeddings <- get_embeddings(loom)

close_loom(loom)

# load motif enrichment results
motifEnrichment <- read.csv(file = "./outs_2023-04-12/pySCENIC_CTX_regulons.csv", header = TRUE, skip = 1)
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

### Add desired data to original Seurat object (RNA only)

seurat_rna <- readRDS(seurat_path)

# Extract underlying matrix from regulonAUC object
AUCmat <- AUCell::getAUC(regulonAUC)
AUCmat[1:10, 1:10]
rownames(AUCmat)

# Add as assay to Seurat object
seurat_rna[['pyscenicAUC']] <- CreateAssayObject(data = AUCmat)
saveRDS(seurat_rna, file = paste0(out_dir, "/", analysis, "_pyscenic2seurat.rds"))

# Extract cell cluster information from Seurat obj
cellClusters <- as.data.frame(Idents(seurat_rna))
head(cellClusters)
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters$`Idents(seurat_rna)`) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
pdf(file = "./outs_2023-04-12/pyscenic_regulon_activity_by_celltype.pdf",
    width = 10,
    height = 12)
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
dev.off()
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later
# To see the exact values
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))

# To identify cluster-specific regulons
# To calculate and visualize the RSS(Regulon Specificity Score)
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters[colnames(regulonAUC),])
## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
fig <- ggplotly(rssPlot$plot)
htmlwidgets::saveWidget(as_widget(fig), 
                        "./outs_2023-04-12/pyscenic_cluster_specific_regulons.html")
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
CellType <- levels(Idents(seurat_rna))
for (i in CellType) {
  pdf(paste0("./outs_2023-04-12/", gsub("/", "_", i), "_specific_regulon.pdf"),
    width = 7,
    height = 7)
  p <- plotRSS_oneSet(rss, setName = i) # cluster ID
  print(p)
  dev.off()
}

# Cell states based on the GRN activity
# cluster the cells based on regulon activity
# SCENIC provides some wrapper functions to get a quick overview and detect potential associations of cell types or states with specific regulons.
# projecting the AUC and TF expression onto 2D projections (such as t-SNE and Umap)
# List of embeddings available:
cat(names(embeddings), sep="\n")
# Overview of these embeddings (see below for details)
regulonsToPlot <- "Sox9_(+)"
options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2, ceiling(length(names(embeddings))/2)))
for (selectedEmbedding in names(embeddings)){
  AUCell::AUCell_plotTSNE(embeddings[[selectedEmbedding]], exprMat, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, 
                          sub=selectedEmbedding)
}
  
