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
so <- readRDS(file.path("F:\\dataset20230202/", "dataset_seu_20PC.rds"))
sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

# Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

# --------------------------------------------------------------
# select resolution = 0.2
so <- SetIdent(so, value = "integrated_snn_res.0.2")
so@meta.data$cluster_id <- Idents(so)
sce$cluster_id <- Idents(so)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "F:\\dataset20230202\\dataset_so_PC20_res0.2_cluster_numbers.csv")

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
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\cell frequency for each cluster.png",
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
cs <- sample(colnames(so), 5e3)
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
  png(paste0("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\", "DR_", id, "_res0.2.png"),
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
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\vlnplot-qc_all.png",
    width = 1200,
    height = 600)
VlnPlot(object = so, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb"), ncol = 4, pt.size = 0)
dev.off()
# vlnplot-featureRNA
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\vlnplot-featureRNA.png",
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
write.csv(feature_by_sample_table, "F:\\dataset20230202\\dataset_so_PC20_res0.2_cluster_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so$nCount_RNA, row.names = so$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "F:\\dataset20230202\\dataset_so_PC20_res0.2_cluster_QC_count_by_sample.csv")

# Generate summary statistics per cluster
char <- as.character(so$cluster_id)
feature_by_cluster <- as.data.frame(so$nFeature_RNA, row.names = char)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "F:\\dataset20230202\\dataset_so_PC20_res0.2_cluster_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so$nCount_RNA, row.names = char)
count_by_cluster_table <- describeBy(count_by_cluster, group = char, mat = TRUE)
write.csv(count_by_cluster_table, "F:\\dataset20230202\\dataset_so_PC20_res0.2_cluster_QC_count_by_cluster.csv")

# DETERMINE HOW MANY CELL TYPES ARE PRESENT IN DATASET BEFORE MAKING NEW SEURAT OBJECTS
library(unikn)
DefaultAssay(so) <- "RNA"
# Dimplot-clusters
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Dimplot-clusters.png",
    width = 1200,
    height = 1200)
DimPlot(so, reduction = "umap", pt.size = 0.001, label = TRUE, label.size = 4) + theme(aspect.ratio = 1) + scale_color_manual(values = usecol("pal_unikn_pair", 40))
dev.off()

# Dimplot-split_by_group
so$group_id <- factor(x = so$group_id, levels = c("ZT6", "ZT18", "SD6h ZT0-ZT6", "SD3D 8h/day (ZT0-ZT8)", "SD3D R7"))
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Dimplot-split_by_group.png",
    width = 4000,
    height = 800)
DimPlot(so, reduction = "umap", split.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1) + scale_color_manual(values = usecol("pal_unikn_pair", 40))
dev.off()

# Find all markers
DefaultAssay(so) <- "RNA"
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers, "F:\\dataset20230202\\dataset_so_PC20_res0.2_cluster_genes-RNA.csv")

# Markers
#Oligo
# Oligo_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Oligo_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Plp1","Mobp","Mog","Cnp","St18","Mag","Sox6","Cldn11", "Olig1","Lhfpl3", "Aspa", "Ermn", "Opalin"))
dev.off()
# Oligo_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Oligo_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Plp1","Mobp","Mog","Cnp","St18","Mag","Sox6","Cldn11", "Olig1","Lhfpl3", "Aspa", "Ermn", "Opalin"), reduction = "umap", coord.fixed = 1)
dev.off()
# Oligo precursor cell
# Oligo_precursor_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Oligo_precursor_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Pdgfra", "Sox6", "Olig1"))
dev.off()
# Oligo_precursor_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Oligo_precursor_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Pdgfra", "Sox6", "Olig1"), reduction = "umap", coord.fixed = 1)
dev.off()
# Oligo progenitor cell
# Oligo_progenitor_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Oligo_progenitor_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Pdgfra", "Vcan"))
dev.off()
# Oligo_progenitor_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Oligo_progenitor_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Pdgfra", "Vcan"), reduction = "umap", coord.fixed = 1)
dev.off()
#ASTRO
# Astro_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Astro_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Slc1a3","Gpr37l1","Gpc5","Mfge8", "Acsbg1", "Aqp4", "Cldn10", "Gja1"))
dev.off()
# Astro_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Astro_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Slc1a3","Gpr37l1","Gpc5","Mfge8", "Acsbg1", "Aqp4", "Cldn10", "Gja1"), reduction = "umap", coord.fixed = 1)
dev.off()
#MICRO
# Micro_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Micro_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Ctss","P2ry12","C1qa","C1qc","C1qb","Tgfbr1", "Cx3cr1"))
dev.off()
# Micro_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Micro_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Ctss","P2ry12","C1qa","C1qc","C1qb","Tgfbr1", "Cx3cr1"), reduction = "umap", coord.fixed = 1)
dev.off()
#Endo
# Endo_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Endo_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Flt1","Cldn5","Ly6a","Igfbp7", "Ly6c1", "Slco1a4"))
dev.off()
# Endo_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Endo_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Flt1","Cldn5","Ly6a","Igfbp7", "Ly6c1", "Slco1a4"), reduction = "umap", coord.fixed = 1)
dev.off()
# Neuron
# Neuron_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Neuron_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Kcnq5", "Map2", "Snhg11", "Gria2", "Meg3", "Ppp2r2b", "Gria3", "Snap25"))
dev.off()
# Neuron_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Neuron_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Kcnq5", "Map2", "Snhg11", "Gria2", "Meg3", "Ppp2r2b", "Gria3", "Snap25"), reduction = "umap", coord.fixed = 1)
dev.off()
# Inh-neuron
# Inh-neuron_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-neuron_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Gad2"))
dev.off()
# Inh-neuron_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-neuron_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad2"), reduction = "umap", coord.fixed = 1)
dev.off()
# Exc-neuron
# Exc-neuron_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-neuron_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Neurod6", "Nm1", "Slc17a7", "Sv2b"))
dev.off()
# Exc-neuron_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-neuron_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Neurod6", "Nm1", "Slc17a7", "Sv2b"), reduction = "umap", coord.fixed = 1)
dev.off()
# Interneuron
# Interneuron_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Interneuron_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Gad1", "Gad2", "Slc32a1"))
dev.off()
# Interneuron_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Interneuron_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad1", "Gad2", "Slc32a1"), reduction = "umap", coord.fixed = 1)
dev.off()
# Inh-PV
# Inh-PV_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-PV_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Erbb4","Pvalb"))
dev.off()
# Inh-PV_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-PV_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Erbb4","Pvalb"), reduction = "umap", coord.fixed = 1)
dev.off()
# Inh-VIP
# Inh-VIP_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-VIP_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Gad2", "Gad1", "Npy", "Vip", "Adarb2"))
dev.off()
# Inh-VIP_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-VIP_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad2", "Gad1", "Npy", "Vip", "Adarb2"), reduction = "umap", coord.fixed = 1)
dev.off()
#Inh-SST
# Inh-SST_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-SST_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Sst"))
dev.off()
# Inh-SST_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-SST_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Sst"), reduction = "umap", coord.fixed = 1)
dev.off()
#Inh-Lamp5
# Inh-Lamp5_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-Lamp5_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Erbb4","Npy","Adarb2","Lamp5","Reln","Slc32a1"))
dev.off()
# Inh-Lamp5_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Inh-Lamp5_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Erbb4","Npy","Adarb2","Lamp5","Reln","Slc32a1"), reduction = "umap", coord.fixed = 1)
dev.off()
# Exc-L5_or_6N
# Exc-L5_or_6N_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L5_or_6N_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Tshz2","Vwc2l"))
dev.off()
# Exc-L5_or_6N_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L5_or_6N_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Tshz2","Vwc2l"), reduction = "umap", coord.fixed = 1)
dev.off()
#Exc-L6CT
# Exc-L6CT_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L6CT_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Foxp2","Hs3st4","Zfpm2"))
dev.off()
# Exc-L6CT_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L6CT_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Foxp2","Hs3st4","Zfpm2"), reduction = "umap", coord.fixed = 1)
dev.off()
# Exc-L2&L4_or_5IT
# Exc-L2&L4_or_5IT_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L2&L4_or_5IT_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Cux2"))
dev.off()
# Exc-L2&L4_or_5IT_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L2&L4_or_5IT_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Cux2"), reduction = "umap", coord.fixed = 1)
dev.off()
# Exc-L2_or_3IT
# Exc-L2_or_3IT_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L2_or_3IT_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Ndst4","Brinp1","Abi3bp"))
dev.off()
# Exc-L2_or_3IT_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L2_or_3IT_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Ndst4","Brinp1","Abi3bp"), reduction = "umap", coord.fixed = 1)
dev.off()
# Exc-L6b(not found)
# Exc-L6b_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L6b_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Lcp1"))
dev.off()
# Exc-L6b_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-L6b_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Lcp1"), reduction = "umap", coord.fixed = 1)
dev.off()
# Exc-Car3
# Exc-Car3_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-Car3_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Col11a1"))
dev.off()
# Exc-Car3_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Exc-Car3_marker_umap_res0.2.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Col11a1"), reduction = "umap", coord.fixed = 1)
dev.off()
# Macrophage
# Macrophage_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Macrophage_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Dab2"))
dev.off()
# Macrophage_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Macrophage_marker_umap_res0.2.png",
    width = 1200,
    height = 600)
FeaturePlot(so, features = c("Dab2"), reduction = "umap", coord.fixed = 1)
dev.off()
# T cell
# T_cell_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\T_cell_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Cd44","Nrp1"))
dev.off()
# T_cell_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\T_cell_marker_umap_res0.2.png",
    width = 1200,
    height = 600)
FeaturePlot(so, features = c("Cd44","Nrp1"), reduction = "umap", coord.fixed = 1)
dev.off()

# Glial cell
# Glial cell_marker_dot_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Glial cell_marker_dot_res0.2.png",
    width = 1200,
    height = 600)
DotPlot(so, features = c("Egfr","Hes1","Plp1","Sox10","Tnc","Mki67"))
dev.off()
# Glial cell_marker_umap_res0.2
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\Glial cell_marker_umap_res0.2.png",
    width = 1200,
    height = 600)
FeaturePlot(so, features = c("Egfr","Hes1","Plp1","Sox10","Tnc","Mki67"), reduction = "umap", coord.fixed = 1)
dev.off()

# Assign cell type identity to clusters
so.renamed <- RenameIdents(so,
                           `0` = "Exc-L2&L4/5IT", `1` = "Astro", `2` = "Exc-neuron1", `3` = "Exc-neuron2",
                           `4` = "Oligo", `5` = "Micro", `6` = "Exc-L6CT", `7` = "Inh-PV&SST", 
                           `8` = "Inh-VIP&Lamp5", `9` = "OPC", `10` = "Neuron1(Map2+)", `11` = "Exc-neuron3",
                           `12` = "Exc-neuron4", `13` = "Neuron2", `14` = "Glial cell", `15` = "Astro(Slc1a3+, Dab2+)",
                           `16` = "Astro(Gja1+)", `17` = "Neuron3", `18` = "Inh-neuron1", `19` = "Endo(Igfbp7+)",
                           `20` = "Endo", `21` = "Exc-Car3")
DefaultAssay(so.renamed) <- "RNA"
# rename_clusters_dimplot
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\rename_clusters_dimplot.png",
    width = 1200,
    height = 1200)
DimPlot(so.renamed, reduction = "umap",label = TRUE, label.size = 10) + theme(aspect.ratio = 1)
dev.off()

# 3.40-rename_clusters_dimplot_check_difference_among_group
# DimPlot(so.renamed, reduction = "umap", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)

# rename_clusters_dimplot_check_difference_among_group_1
png("E:\\ZZY\\snRNAseq_analysis\\20230202-analyze_with_new_code\\#3-Annotation\\3.2#-Annotation-res0.2\\rename_clusters_dimplot_check_difference_among_group_1.png",
    width = 4000,
    height = 800)
DimPlot(so.renamed, reduction = "umap", split.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
dev.off()

saveRDS(so.renamed, file.path("F:\\dataset20230202/", "so.renamed_res0.2.rds"))


