# slingshot procedure: http://barcwiki.wi.mit.edu/wiki/SOP/scRNA-seq/Slingshot
# gam procedure: https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html#slingshot
# loading packages
library(Seurat) 
library(SingleCellExperiment)
library(slingshot)

# loading data
seu <- readRDS("../../#9-subset Inh neuron/so_inh.renamed_res0.1.rds")
groups <- c("ZT6", "ZT18","SD3D 8h/day (ZT0-ZT8)")
seu <- subset(x = seu, idents = c("Inh-VIP(Penk+)","Inh-Lamp5(Dbp+)","Inh-Adarb2+Lhx6-Pbx3+","Inh-PAX6"), subset = group_id %in% groups, invert = FALSE)
# cells_by_sample <- split(colnames(seu), seu$sample_id)
# seu <- lapply(cells_by_sample, function(i) {
#   subset(seu, cells = i)
# })

# analyze NS,NW,SD3D
sce <- as.SingleCellExperiment(seu, assay = "RNA")
#this has the cell classification
table(sce$ident)

sce <- slingshot(sce, clusterLabels = sce$ident, reducedDim = "UMAP",
                 allow.breaks = FALSE)
summary(sce$slingPseudotime_1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   23.97   31.77   35.28   48.25   87.26 

# get the lineages:
lnes <- getLineages(reducedDim(sce,"UMAP"), sce$ident)
lnes@metadata$lineages

# visualize the pseudotime and lineage
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(dplyr)
# this define the cluster color. You can change it with different color scheme.
my_color <- createPalette(length(levels(sce$ident)), c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(sce$ident))
count <- colData(sce)
slingshot_df <- data.frame(slingPseudotime_1=count$slingPseudotime_1, ident=count$ident)

# re-order y-axis for better figure: This should be tailored with your own cluster names
slingshot_df$ident = factor(slingshot_df$ident, levels=c("Inh-Lamp5(Dbp+)","Inh-VIP(Penk+)","Inh-PAX6","Inh-Adarb2+Lhx6-Pbx3+"))

pdf("./slingshot_neurons_pseudotime_N1_N2_N4_adarb2.pdf")
ggplot(slingshot_df, aes(x = slingPseudotime_1, y = ident, 
                         colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)
dev.off()

pdf("./slingshot_neurons_N1_N2_N4_adarb2.pdf",
    width = 10,
    height = 10,
    onefile = F)
plot(reducedDims(sce)$UMAP, col = my_color[as.character(sce$ident)], 
     pch=16, 
     asp = 1)
legend("bottomleft",legend = names(my_color[levels(sce$ident)]),  
       fill = my_color[levels(sce$ident)])
lines(SlingshotDataSet(lnes), lwd=2, type = 'lineages', col = c("black"))
dev.off()

# GAM general additive model for identifying temporally expressed genes
library(gam)
t <- sce$slingPseudotime_1
# for time, only look at the 100 most variable genes 
Y <- log1p(assay(sce,"logcounts"))

var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- gam(z ~ lo(t), data=d)
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

## Plot the top 100 genes' expression
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]

heatdata <- assays(sce)$logcounts[topgenes, order(t, na.last = NA)]
heatdata <- as.matrix(heatdata)
heatclus <- sce$ident[order(t, na.last = NA)]

pdf("./slingshot_neurons_N1_N2_N4_adarb2_gene_heatmap.pdf",
    width = 10,
    height = 20,
    onefile = F)
heatmap(heatdata, Colv = NA,
        ColSideColors = my_color[heatclus],cexRow = 0.5,cexCol = 1,
        labCol = NA)
dev.off()

# check branch genes
seu <- readRDS("../../#9-subset Inh neuron/so_inh.renamed_res0.1.rds")
png("./Kcnb2_marker_dot_res0.1.png",
    width = 1200,
    height = 600)
DotPlot(seu, features = c("Kcnb2", "Vip", "Prkg1", "Pax6", "Pbx3", "Meis2", "Pcbp3", "Dlx6os1"))
dev.off()
# Adarb2_Lhx6_marker_umap_res0.1
png("./Kcnb2_marker_umap_res0.1.png",
    width = 1200,
    height = 1200)
FeaturePlot(seu, features = c("Kcnb2", "Vip", "Prkg1", "Pax6" ,"Pbx3", "Meis2", "Pcbp3", "Dlx6os1"), reduction = "umap", coord.fixed = 1)
dev.off()
png("./marker_umap_res0.1.png",
    width = 1200,
    height = 1200)
FeaturePlot(seu, features = c("Adarb2","Lhx6","Lamp5","Pax6","Pvalb","Sst","Vip","Pbx3"), reduction = "umap", coord.fixed = 1)
dev.off()