# ---------------------------------------------------
## DR
## loading packages
library(so_inhrat)
library(cowplot)
library(SingleCellExperiment)
library(ggplot2)

## loading data
so.renamed <- readRDS("so.renamed_res0.2.rds")
so_inh <- subset(x = so.renamed, idents = c("Inh-PV&SST"), invert = FALSE)

## reclustering
#Convert SCE to so_inhratObject
sce <- as.SingleCellExperiment(so_inh, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_inh <- lapply(cells_by_sample, function(i){
  subset(so_inh, cells = i)
})

#Normalize, find variable genes, and scale
so_inh <- lapply(so_inh, NormalizeData, verbose = FALSE)
so_inh <- lapply(so_inh, FindVariableFeatures, nfeatures = 2e3,
                   selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_inh <- lapply(so_inh, ScaleData, verbose = FALSE)

#Find anchors and integrate 
##Decrease k.filter to minimize number of astrocytes identified per sample
as <- FindIntegrationAnchors(so_inh, verbose = TRUE, k.filter = 125)
so_inh <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_inh) <- "integrated"
so_inh <- ScaleData(so_inh, display.progress = FALSE)

#DIMENSION REDUCTION
so_inh <- RunPCA(so_inh, features = VariableFeatures(object = so_inh), verbose = FALSE)
print(so_inh[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(so_inh,dims = 1:2,reduction = "pca")
DimPlot(so_inh,reduction = "pca")
DimHeatmap(so_inh,dims = 1,cells = 500,balanced = TRUE)
# DimHeatmap: check main so_inhurce of heterogenity in dataset
png("./DimHeatmap.png",
    width = 1800,
    height = 4000)
DimHeatmap(so_inh,dims = 1:30,cells=500,balanced = TRUE)
dev.off()
# determine the PC of dataset
# Method 1: JackStraw
so_inh <- JackStraw(so_inh, num.replicate = 100)
so_inh <- ScoreJackStraw(so_inh, dims = 1:20)
saveRDS(so_inh, file = "./so_inh_temp.RDS")
# Visualization
png("./JackStrawPlot.png",
    width = 1200,
    height = 1200)
JackStrawPlot(so_inh, dims = 1:20)
dev.off()
# Method 2: Elbow
png("./ElbowPlot.png",
    width = 1200,
    height = 1200)
ElbowPlot(so_inh, ndims = 50)
dev.off()
# Calculate the inflection point
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# Determine percent of variation associated with each PC
pct <- so_inh[["pca"]]@stdev / sum(so_inh[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
##43

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
##19

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs
##19

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
# Choosing PC = 19

#Update number of PCs used
so_inh <- RunTSNE(so_inh, reduction = "pca", dims = seq_len(19),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_inh <- RunUMAP(so_inh, reduction = "pca", dims = seq_len(19),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_inh <- FindNeighbors(so_inh, reduction = "pca", dims = seq_len(19), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_inh <- FindClusters(so_inh, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_inh, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_inh, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
png("./DR_inh.png",
    width = 1200,
    height = 1200)
plot_grid(plotlist = ps, ncol = 1)
dev.off()

#Save SeuratObject
saveRDS(so_inh, file.path("so_inh_19PC.rds"))

# ------------------------------------------------------------
## Annotation
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

# Load data and convert to SCE
so <- readRDS(file.path("./", "so_inh_19PC.rds"))
sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

# Define resolution
# Check number of clusters in each resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

# select resolution = 0.2
so <- SetIdent(so, value = "integrated_snn_res.0.2")
so@meta.data$cluster_id <- Idents(so)
sce$cluster_id <- Idents(so)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "./dataset_so_PC19_res0.2_cluster_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

# Choose color palettes for cluster, sample, group IDs, and # cells
library(CATALYST)
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green", "purple"), gids)

# calculate frequency of cell numbers in each cluster for each sample
# sum up all cell numbers in each sample, and cell number in each cluster is divided by the sum
fqs <- prop.table(n_cells, margin = 2)
mat <- as.matrix(unclass(fqs))
# cell frequency for each cluster
png("./cell frequency for each cluster.png",
    width = 600,
    height = 800)
Heatmap(mat,
        col = rev(brewer.pal(11, "RdGy")[-6]),
        name = "Frequency",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_title = "cluster_id",
        column_title = "sample_id",
        column_title_side = "bottom",
        rect_gp = gpar(col = "white"),
        cell_fun = function(i, j, x, y, width, height, fill) {
          grid.text(round(mat[j, i] * 100, 2),
                    x = x, y = y,
                    gp = gpar(col = "white", fontsize = 8)
          )
        }
)
dev.off()

# DR
cs <- sample(colnames(so), 2e3)
.plot_dr <- function(so, dr, id) {
  DimPlot(so, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
    scale_color_manual(id, values = get(paste0(id, "_pal"))) +
    guides(col = guide_legend(
      nrow = 10,
      override.aes = list(size = 3, alpha = 1)
    )) +
    theme_void() + theme(aspect.ratio = 1)
}

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  png(paste0("./", "DR_", id, "_res0.2.png"),
      width = 600,
      height = 600)
  cat("## ", id, "\n")
  p1 <- .plot_dr(so, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
  dev.off()
}
# DR_cluster_id_res0.2
# DR_group_id_res0.2
# DR_sample_id_res0.2

# QC METRICS CHECK
mito.genes <- grep(pattern = "^mt-", x = rownames(so@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so@assays[["RNA"]][mito.genes, ]) / Matrix::colSums(so@assays[["RNA"]])
so$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^Rp[sl]", x = rownames(so@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so@assays[["RNA"]][rb.genes, ]) / Matrix::colSums(so@assays[["RNA"]])
so$percent.rb <- percent.rb

# vlnplot-qc_all
png("./vlnplot-qc_all.png",
    width = 1200,
    height = 600)
VlnPlot(object = so, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb"), ncol = 4, pt.size = 0)
dev.off()
# vlnplot-featureRNA
png("./vlnplot-featureRNA.png",
    width = 1200,
    height = 600)
VlnPlot(object = so, features = c("nFeature_RNA"), pt.size = 0) + stat_summary(fun.y = median, geom = "point", shape = 23, size = 2)
dev.off()

# Generate summary statistics for entire dataset
summary(so$nFeature_RNA)
summary(so$nCount_RNA)

library(pastecs)
stat.desc(so$nFeature_RNA)
stat.desc(so$nCount_RNA)

library(psych)
describe(so$nFeature_RNA)
describe(so$nCount_RNA)

# Generate summary statistics per sample
library(data.table)
library(psych)
feature_by_sample <- as.data.frame(so$nFeature_RNA, row.names = so$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "./dataset_so_PC19_res0.2_cluster_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so$nCount_RNA, row.names = so$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "./dataset_so_PC19_res0.2_cluster_QC_count_by_sample.csv")

# Generate summary statistics per cluster
char <- as.character(so$cluster_id)
feature_by_cluster <- as.data.frame(so$nFeature_RNA, row.names = char)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "./dataset_so_PC19_res0.2_cluster_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so$nCount_RNA, row.names = char)
count_by_cluster_table <- describeBy(count_by_cluster, group = char, mat = TRUE)
write.csv(count_by_cluster_table, "./dataset_so_PC19_res0.2_cluster_QC_count_by_cluster.csv")

# DETERMINE HOW MANY CELL TYPES ARE PRESENT IN DATASET BEFORE MAKING NEW SEURAT OBJECTS
library(unikn)
DefaultAssay(so) <- "RNA"
# Dimplot-clusters
png("./Dimplot-clusters.png",
    width = 1200,
    height = 1200)
DimPlot(so, reduction = "umap", pt.size = 0.001, label = TRUE, label.size = 4) + theme(aspect.ratio = 1) + scale_color_manual(values = usecol("pal_unikn_pair", 40))
dev.off()

# Dimplot-split_by_group
so$group_id <- factor(x = so$group_id, levels = c("ZT6", "ZT18", "SD6h ZT0-ZT6", "SD3D 8h/day (ZT0-ZT8)", "SD3D R7"))
png("./Dimplot-split_by_group.png",
    width = 4000,
    height = 800)
DimPlot(so, reduction = "umap", split.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1) + scale_color_manual(values = usecol("pal_unikn_pair", 40))
dev.off()

# Find all markers
DefaultAssay(so) <- "RNA"
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers, "./dataset_so_PC19_res0.2_cluster_genes-RNA.csv")

# Markers
# Inh-neuron
# Inh-neuron_marker_dot_res0.2
png("./Inh-neuron_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Gad2"))
dev.off()
# Inh-neuron_marker_umap_res0.2
png("./Inh-neuron_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad2"), reduction = "umap", coord.fixed = 1)
dev.off()
# Interneuron
# Interneuron_marker_dot_res0.2
png("./Interneuron_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Gad1", "Gad2", "Slc32a1"))
dev.off()
# Interneuron_marker_umap_res0.2
png("./Interneuron_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad1", "Gad2", "Slc32a1"), reduction = "umap", coord.fixed = 1)
dev.off()
# Inh-PV
# Inh-PV_marker_dot_res0.2
png("./Inh-PV_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Erbb4","Pvalb"))
dev.off()
# Inh-PV_marker_umap_res0.2
png("./Inh-PV_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Erbb4","Pvalb"), reduction = "umap", coord.fixed = 1)
dev.off()
# Inh-VIP
# Inh-VIP_marker_dot_res0.2
png("./Inh-VIP_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Gad2", "Gad1", "Npy", "Vip", "Adarb2"))
dev.off()
# Inh-VIP_marker_umap_res0.2
png("./Inh-VIP_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad2", "Gad1", "Npy", "Vip", "Adarb2"), reduction = "umap", coord.fixed = 1)
dev.off()
#Inh-SST
# Inh-SST_marker_dot_res0.2
png("./Inh-SST_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Sst"))
dev.off()
# Inh-SST_marker_umap_res0.2
png("./Inh-SST_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Sst"), reduction = "umap", coord.fixed = 1)
dev.off()

# Assign cell type identity to clusters
so.renamed <- RenameIdents(so,
                           `0` = "Inh-PV", `1` = "Inh-SST", `2` = "Inh-PV", `3` = "Inh-SST",
                           `4` = "Inh-SST", `5` = "Inh-PV", `6` = "else", `7` = "Inh-SST", 
                           `8` = "else", `9` = "else")
DefaultAssay(so.renamed) <- "RNA"
# rename_clusters_dimplot
png("./rename_clusters_dimplot.png",
    width = 1200,
    height = 1200)
DimPlot(so.renamed, reduction = "umap",label = TRUE, label.size = 10) + theme(aspect.ratio = 1)
dev.off()

# 3.40-rename_clusters_dimplot_check_difference_among_group
# DimPlot(so.renamed, reduction = "umap", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)

# rename_clusters_dimplot_check_difference_among_group_1
png("./rename_clusters_dimplot_check_difference_among_group_1.png",
    width = 4000,
    height = 800)
DimPlot(so.renamed, reduction = "umap", split.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
dev.off()

saveRDS(so.renamed, file.path("./", "so_inh.renamed_res0.2.rds"))
