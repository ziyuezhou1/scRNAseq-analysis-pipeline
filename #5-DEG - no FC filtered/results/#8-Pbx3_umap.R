# loading packages
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# loading data
so <- readRDS("../data/so.renamed_res0.2.rds")

png("./Pbx3_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Pbx3"), reduction = "umap", coord.fixed = 1) +
  # scale_color_continuous(colours = )
        theme(text = element_text(size = 40),
        axis.text = element_text(size = 30),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()

png("./PV_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Pvalb"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()

png("./Sst_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Sst"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()

png("./Lamp5_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Lamp5"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()
png("./Vip_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Vip"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()
png("./Gad1_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad1"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()
png("./Gad1_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Gad1"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()
png("./Meis2_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Meis2"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()
png("./Pcbp3_umap.png",
    width = 1200,
    height = 1200)
FeaturePlot(so, features = c("Pcbp3"), reduction = "umap", coord.fixed = 1) +
  theme(text = element_text(size = 50),
        axis.text = element_text(size = 40),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(2, 'cm'))
dev.off()