# Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

# Load seu
seu <- readRDS(file.path("../#1-QC/", "seu.rds"))

# Split by sample
# cells_by_sample <- split(colnames(seu), seu$sample_id)
# seu <- lapply(cells_by_sample, function(i) {
#   subset(seu, cells = i)
# })

# Normalize, find variable genes, and scale
# seu <- lapply(seu, NormalizeData, verbose = FALSE)
# seu <- lapply(seu, FindVariableFeatures,
#              nfeatures = 2000,
#              selection.method = "vst", do.plot = FALSE, verbose = FALSE
# )
# seu <- lapply(seu, ScaleData, verbose = FALSE)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000,
                            selection.method = "vst",
                            do.plot = FALSE,
                            verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)

# DIMENSION REDUCTION
seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = FALSE)
print(seu[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seu,dims = 1:2,reduction = "pca")
DimPlot(seu,reduction = "pca")
DimHeatmap(seu,dims = 1,cells = 500,balanced = TRUE)
# DimHeatmap: check main seuurce of heterogenity in dataset
png("./DimHeatmap.png",
    width = 1800,
    height = 4000)
DimHeatmap(seu,dims = 1:30,cells=500,balanced = TRUE)
dev.off()
# determine the PC of dataset
# Method 1: JackStraw
seu <- JackStraw(seu, num.replicate = 100)
seu <- ScoreJackStraw(seu, dims = 1:20)
saveRDS(seu, file = "F:\\dataset20230202\\seu_temp.RDS")
# Visualization
png("./JackStrawPlot.png",
    width = 1200,
    height = 1200)
JackStrawPlot(seu, dims = 1:20)
dev.off()
# Method 2: Elbow
png("./ElbowPlot.png",
    width = 1200,
    height = 1200)
ElbowPlot(seu, ndims = 50)
dev.off()
# Calculate the inflection point
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# Determine percent of variation associated with each PC
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
##41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
##20

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs
##20

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize
png("./ElbowPlot_with_values.png",
    width = 1200,
    height = 1200)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
# Choosing PC = 20

# Update number of PCs used
seu <- RunTSNE(seu,
              reduction = "pca", dims = seq_len(20),
              seed.use = 1, do.fast = TRUE, verbose = FALSE
)
seu <- RunUMAP(seu,
              reduction = "pca", dims = seq_len(20),
              seed.use = 1, verbose = FALSE
)

# CLUSTERING
seu <- FindNeighbors(seu, reduction = "pca", dims = seq_len(20), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) {
  seu <- FindClusters(seu, resolution = res, random.seed = 1, verbose = FALSE)
}

# DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(seu, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(seu, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd,
            nrow = 1,
            rel_widths = c(1, 1, 0.5)
  )
})
# 2.3-DR
png("./DR.png",
    width = 1200,
    height = 1200)
plot_grid(plotlist = ps, ncol = 1)
dev.off()

# SAVE SeuratObject
saveRDS(seu, file.path("./", "seu_20PC.rds"))
