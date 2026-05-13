# Single-cell RNA-seq analysis
# Load libraries
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(tidyr)
library(ggplot2)
library(Seurat)

# Read in the SingleCellExperiment object with the filtered raw counts
seurat <- readRDS(".\\data\\so.renamed_res0.2.rds")
# ------------------------------------------------------------------
# N5 vs N4
# Create a list object to store results
if (!file.exists("results/N5_vs_N4")) dir.create("results/N5_vs_N4")
N5_VS_N4 <- list()

for (i in levels(Idents(seurat))){
  N5_VS_N4[[i]] <- FindMarkers(object = seurat,
                               group.by = "group_id",
                               ident.1 = "SD3D R7",
                               ident.2 = "SD3D 8h/day (ZT0-ZT8)",
                               subset.ident = i,
                               test.use = "wilcox")
}
head(N5_VS_N4)
to_name_csv <- levels(Idents(seurat))
to_name_csv <- gsub("/", "or", to_name_csv)
for (i in seq.int(1,23,1)){
  write.csv(N5_VS_N4[[i]], file = paste0(".\\results\\N5_vs_N4\\N5_vs_N4_", to_name_csv[[i]], ".csv"))
}
df <- as.data.frame(do.call(rbind, N5_VS_N4))
df$names <- row.names(df)
df <- df %>% 
  separate(col = names, into = c("groups", "Genes"), sep = "\\.")
df$threshold <- "non-sig" # p > 0.05
df[(df$avg_log2FC > 0) & (df$p_val_adj <= 0.05),]$threshold <- "sig log2FC > 0" # sigFC > 1
df[(df$avg_log2FC < 0) & (df$p_val_adj <= 0.05),]$threshold <- "sig log2FC < 0" # sigFC < 1

png("./figures/N5_vs_N4.png",
    width = 2500,
    height = 800)
p <- ggplot(df, aes(x = groups, y = avg_log2FC, color = threshold)) +
  geom_jitter(size = 2, width = 0.2)
p + scale_color_manual(values=c("#999999", "#386087", "#9A1200")) +
  theme_classic() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  labs(title = "DEG N5 vs. N4", x = "Clusters", y = "log2FC(N5/N4)")
dev.off()

# merge subtypes into general types
Idents(seurat) <- seurat$RNA_snn_res.0.2
seurat.renamed <- RenameIdents(seurat,
                               `0` = "Exc-neuron", `1` = "Astro", `2` = "Oligo", `3` = "Exc-neuron",
                               `4` = "Micro", `5` = "Exc-neuron", `6` = "Inh-neuron", `7` = "Exc-neuron", 
                               `8` = "Inh-neuron", `9` = "OPC", `10` = "Exc-neuron", `11` = "Exc-neuron",
                               `12` = "Exc-neuron", `13` = "Inh-neuron", `14` = "Astro", `15` = "Exc-neuron",
                               `16` = "Glial cell(Egfr+)", `17` = "Endo", `18` = "Endo", `19` = "Astro",
                               `20` = "Exc-neuron", `21` = "Oligo", `22` = "OPC")
levels(Idents(seurat.renamed))
N5_VS_N4_re <- list()

for (i in levels(Idents(seurat.renamed))){
  N5_VS_N4_re[[i]] <- FindMarkers(object = seurat.renamed,
                                  group.by = "group_id",
                                  ident.1 = "SD3D R7",
                                  ident.2 = "SD3D 8h/day (ZT0-ZT8)",
                                  subset.ident = i,
                                  test.use = "wilcox")
}
head(N5_VS_N4_re)
to_name_csv <- levels(Idents(seurat.renamed))
if (!file.exists("results/N5_vs_N4_re")) dir.create("results/N5_vs_N4_re")
for (i in seq.int(1,8,1)){
  write.csv(N5_VS_N4_re[[i]], file = paste0(".\\results\\N5_vs_N4_re\\N5_vs_N4_re_", to_name_csv[[i]], ".csv"))
}
df_re <- as.data.frame(do.call(rbind, N5_VS_N4_re))
df_re$names <- row.names(df_re)
df_re <- df_re %>% 
  separate(col = names, into = c("groups", "Genes"), sep = "\\.")
df_re$threshold <- "non-sig" # p > 0.05
df_re[(df_re$avg_log2FC > 0.25) & (df_re$p_val_adj <= 0.05),]$threshold <- "sig log2FC > 0.25" # sigFC > 1
df_re[(df_re$avg_log2FC < -0.25) & (df_re$p_val_adj <= 0.05),]$threshold <- "sig log2FC < -0.25" # sigFC < 1
df_re$groups <- factor(df_re$groups, levels = unique(df_re$groups))
df_re$groups <- factor(df_re$groups, levels = c("Exc-neuron", "Inh-neuron", "Astro", "Micro", "Oligo", "OPC", "Glial cell(Egfr+)", "Endo"))
png("./figures/N5_vs_N4_re.png",
    width = 3000,
    height = 600)
p <- ggplot(df_re, aes(x = groups, y = avg_log2FC, color = threshold)) +
  geom_jitter(size = 2, width = 0.2) +
  scale_y_continuous(limits = c(-3, 3))
p + scale_color_manual(values=c("#999999", "#386087", "#9A1200")) +
  theme_classic() +
  theme(text = element_text(size = 50), plot.title = element_text(hjust = 0.5),
        axis.title.x.bottom = element_text(size = 50),
        axis.text.x.bottom = element_text(size = 50)) +
  labs(title = "DEG R7 vs. SD3D", x = "Clusters", y = "log2FC") +
  guides(colour = guide_legend(override.aes = list(size=20)))
dev.off()
