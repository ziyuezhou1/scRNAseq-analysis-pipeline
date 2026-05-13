# loading library
library(dplyr)
library(Seurat)
library(DElegate)
library(ggplot2)
library(ggrepel)

# loading data
seurat <- readRDS("./data/so.renamed_res0.2.rds")

# ZT18 vs. ZT6
de_results <- findDE(object = seurat, group_column = "group_id",
                     compare = c('ZT18', 'ZT6'),
                     method = 'deseq')
de_results$category <- "non-sig"
de_results[de_results$log_fc > 1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc > 1"
de_results[de_results$log_fc < -1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc < -1"

png("./ZT18 vs ZT6 pseudobulk.png",
    width = 1000,
    height = 1000)
ggplot(de_results, aes(x = log_fc, y = -log10(padj), color = category)) +
  geom_point(shape = 16) +
  scale_color_manual(values=c("#999999", "#386087", "#9A1200")) +
  theme(text = element_text(size = 40)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ggtitle('ZT18 vs ZT6')
dev.off()

# SD1D vs. ZT6
options(future.globals.maxSize=1471026299)
de_results <- findDE(object = seurat, group_column = "group_id",
                     compare = c('SD6h ZT0-ZT6', 'ZT6'),
                     method = 'deseq')
write.csv(de_results, file = "./SD1D vs ZT6 pseudobulk.csv")

de_results$category <- "non-sig"
de_results[de_results$log_fc > 1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc > 1"
de_results[de_results$log_fc < -1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc < -1"

png("./SD1D vs ZT6 pseudobulk.png",
    width = 1000,
    height = 1000)
ggplot(de_results, aes(x = log_fc, y = -log10(padj), color = category)) +
  geom_point(shape = 16) +
  scale_color_manual(values=c("#999999", "#386087", "#9A1200")) +
  theme(text = element_text(size = 40)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ggtitle('SD1D vs ZT6')
dev.off()

# SD3D vs. ZT6
de_results <- findDE(object = seurat, group_column = "group_id",
                     compare = c('SD3D 8h/day (ZT0-ZT8)', 'ZT6'),
                     method = 'deseq')
write.csv(de_results, file = "./SD3D vs ZT6 pseudobulk.csv")
de_results$category <- "non-sig"
de_results[de_results$log_fc > 1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc > 1"
de_results[de_results$log_fc < -1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc < -1"

png("./SD3D vs ZT6 pseudobulk.png",
    width = 1000,
    height = 1000)
ggplot(de_results, aes(x = log_fc, y = -log10(padj), color = category)) +
  geom_point(shape = 16) +
  scale_color_manual(values=c("#999999", "#386087", "#9A1200")) +
  theme(text = element_text(size = 40)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ggtitle('SD3D 8h/day (ZT0-ZT8) vs ZT6')
dev.off()

# R7 vs. ZT6
de_results <- findDE(object = seurat, group_column = "group_id",
                     compare = c('SD3D R7', 'ZT6'),
                     method = 'deseq')
write.csv(de_results, file = "./SD3D R7 vs ZT6 pseudobulk.csv")
de_results$category <- "non-sig"
de_results[de_results$log_fc > 1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc > 1"
de_results[de_results$log_fc < -1 & is.na(de_results$padj) == FALSE & de_results$padj < 0.05,]$category <- "sig log_fc < -1"

png("./SD3D R7 vs ZT6 pseudobulk.png",
    width = 1000,
    height = 1000)
ggplot(de_results, aes(x = log_fc, y = -log10(padj), color = category)) +
  geom_point(shape = 16) +
  scale_color_manual(values=c("#999999", "#386087", "#9A1200")) +
  theme(text = element_text(size = 40)) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  ggtitle('SD3D R7 vs ZT6')
dev.off()