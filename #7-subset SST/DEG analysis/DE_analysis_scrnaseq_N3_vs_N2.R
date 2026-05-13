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
seurat <- readRDS(".\\data\\so_inh.renamed_res0.2.rds")
# ------------------------------------------------------------------
# N3 vs N2
# Create a list object to store results
if (!file.exists("results/N3_vs_N2")) dir.create("results/N3_vs_N2")
N3_VS_N2 <- list()

for (i in levels(Idents(seurat))){
  try(
    N3_VS_N2[[i]] <- FindMarkers(object = seurat,
                                 group.by = "group_id",
                                 ident.1 = "SD6h ZT0-ZT6",
                                 ident.2 = "ZT18",
                                 subset.ident = i,
                                 logfc.threshold = 0,
                                 min.pct = -Inf,
                                 min.cells.feature = 0,
                                 min.cells.group = 0,
                                 test.use = "wilcox")
  )
}
head(N3_VS_N2)
to_name_csv <- levels(Idents(seurat))
to_name_csv <- gsub("/", "or", to_name_csv)
for (i in seq.int(1,3,1)){
  write.csv(N3_VS_N2[[i]], file = paste0(".\\results\\N3_vs_N2\\N3_vs_N2_", to_name_csv[[i]], ".csv"))
}
df <- as.data.frame(do.call(rbind, N3_VS_N2))
df$names <- row.names(df)
df <- df %>% 
      separate(col = names, into = c("groups", "Genes"), sep = "\\.")
df$threshold <- "non-sig" # p > 0.05
df[(df$avg_log2FC > 0) & (df$p_val_adj <= 0.05),]$threshold <- "sig log2FC > 0" # sigFC > 1
df[(df$avg_log2FC < 0) & (df$p_val_adj <= 0.05),]$threshold <- "sig log2FC < 0" # sigFC < 1

png("./figures/N3_vs_N2.png",
    width = 2500,
    height = 800)
  p <- ggplot(df, aes(x = groups, y = avg_log2FC, color = threshold)) +
    geom_jitter(size = 2, width = 0.2)
  p + scale_color_manual(values=c("#999999", "#386087", "#9A1200")) +
    theme_classic() +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
    labs(title = "DEG N3 vs. N2", x = "Clusters", y = "log2FC(N3/N2)")
dev.off()
