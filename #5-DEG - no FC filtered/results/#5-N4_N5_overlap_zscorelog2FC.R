# loading package
library(ggplot2)
library(ggrepel)
library(tidyverse)

# create result directory
if (!dir.exists("./N4_N5_overlap_log2FC")) {dir.create("./N4_N5_overlap_log2FC")}

draw_scatter <- function(x){
  # loading data
  N4_vs_N1 <- read.csv(paste0("./N4_vs_N1_re/N4_vs_N1_re_", x, ".csv"))
  N5_vs_N4 <- read.csv(paste0("./N5_vs_N4_re/N5_vs_N4_re_", x, ".csv"))
  # subset padj <= 0.05: only subset out SD3D vs. NS DEG
  N4_vs_N1_sig <- N4_vs_N1[N4_vs_N1$p_val_adj <= 0.05,]
  # changing name
  N4_vs_N1_sig$avg_log2FC_N4_N1_DEGs <- N4_vs_N1_sig$avg_log2FC
  N5_vs_N4$avg_log2FC_N5_N4 <- N5_vs_N4$avg_log2FC
  # find overlap genes: possible rescued DEGs
  overlap <- as.vector(intersect(N4_vs_N1_sig$X, N5_vs_N4$X))
  # subset log2FC out
  N4_vs_N1_overlap <- N4_vs_N1_sig[N4_vs_N1_sig$X %in% overlap,]
  N5_vs_N4_overlap <- N5_vs_N4[N5_vs_N4$X %in% overlap,]
  # merge together
  overlap_df <- merge(N4_vs_N1_overlap,
                      N5_vs_N4_overlap,
                      by = "X")
  overlap_df$N5_N4_sig <- "non-sig"
  overlap_df$color.code <- "#00AFBB"
  if (nrow(overlap_df[overlap_df$p_val_adj.y <= 0.05,]) != 0) {
    overlap_df[overlap_df$p_val_adj.y <= 0.05,]$N5_N4_sig <- "sig"
    overlap_df[overlap_df$N5_N4_sig == "sig",]$color.code <- "#E7B800"
  }
  write.csv(overlap_df, file = paste0("./N4_N5_overlap_log2FC/N4_N5_DEG_overlap_", x, "_log2FC.csv"))
  # draw scatterplot
  png(filename = paste0("./N4_N5_overlap_log2FC/N4_N5_DEG_overlap_", x, "_log2FC.png"),
      width = 1200,
      height = 1200)
  p <- ggplot(overlap_df, aes(x = avg_log2FC_N4_N1_DEGs, y = avg_log2FC_N5_N4,
                              label = X)) +
    annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= "#72be7a", alpha = 0.5) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= "#72be7a", alpha = 0.5) +
    geom_point(aes(colour = N5_N4_sig), size = 5) +
    geom_text_repel(hjust=0, vjust=0, size = 10) +
    theme(text = element_text(size = 50),
          legend.title = element_blank(),
          legend.text = element_text(size = 50)) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_color_manual(values = setNames(overlap_df$color.code, overlap_df$N5_N4_sig)) +
    labs(y = "log2FC(R7/SD3D)", x = "log2FC(SD3D/NS)") +
    guides(color = guide_legend(override.aes = list(size = 20)))
  print(p + ggtitle(x))
  dev.off()
}

celltype <- c("Astro",
              "Endo",
              "Exc-neuron",
              "Glial cell(Egfr+)",
              # "Inh-neuron",
              "Micro",
              "Oligo",
              "OPC")
for (i in celltype) {
  draw_scatter(i)
}

draw_scatter_merge_glial <- function(x){
  # loading data
  N4_vs_N1 <- read.csv(paste0("./N4_vs_N1_re/N4_vs_N1_re_merge_glial_", x, ".csv"))
  N5_vs_N4 <- read.csv(paste0("./N5_vs_N4_re/N5_vs_N4_re_merge_glial_", x, ".csv"))
  # subset padj <= 0.05: only subset out SD3D vs. NS DEG
  N4_vs_N1_sig <- N4_vs_N1[N4_vs_N1$p_val_adj <= 0.05,]
  # changing name
  N4_vs_N1_sig$avg_log2FC_N4_N1_DEGs <- N4_vs_N1_sig$avg_log2FC
  N5_vs_N4$avg_log2FC_N5_N4 <- N5_vs_N4$avg_log2FC
  # find overlap genes: possible rescued DEGs
  overlap <- as.vector(intersect(N4_vs_N1_sig$X, N5_vs_N4$X))
  # subset log2FC out
  N4_vs_N1_overlap <- N4_vs_N1_sig[N4_vs_N1_sig$X %in% overlap,]
  N5_vs_N4_overlap <- N5_vs_N4[N5_vs_N4$X %in% overlap,]
  # merge together
  overlap_df <- merge(N4_vs_N1_overlap,
                      N5_vs_N4_overlap,
                      by = "X")
  overlap_df$N5_N4_sig <- "non-sig"
  overlap_df$color.code <- "#00AFBB"
  if (nrow(overlap_df[overlap_df$p_val_adj.y <= 0.05,]) != 0) {
    overlap_df[overlap_df$p_val_adj.y <= 0.05,]$N5_N4_sig <- "sig"
    overlap_df[overlap_df$N5_N4_sig == "sig",]$color.code <- "#E7B800"
  }
  write.csv(overlap_df, file = paste0("./N4_N5_overlap_log2FC/N4_N5_DEG_overlap_merge_glial_", x, "_log2FC.csv"))
  # draw scatterplot
  png(filename = paste0("./N4_N5_overlap_log2FC/N4_N5_DEG_overlap__merge_glial_", x, "_log2FC.png"),
      width = 1200,
      height = 1200)
  p <- ggplot(overlap_df, aes(x = avg_log2FC_N4_N1_DEGs, y = avg_log2FC_N5_N4,
                              label = X)) +
    annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= "#72be7a", alpha = 0.5) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= "#72be7a", alpha = 0.5) +
    geom_point(aes(colour = N5_N4_sig), size = 5) +
    geom_text_repel(hjust=0, vjust=0, size = 10) +
    theme(text = element_text(size = 50),
          legend.title = element_blank(),
          legend.text = element_text(size = 50)) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_color_manual(values = setNames(overlap_df$color.code, overlap_df$N5_N4_sig)) +
    labs(y = "log2FC(R7/SD3D)", x = "log2FC(SD3D/NS)") +
    guides(color = guide_legend(override.aes = list(size = 20)))
  print(p + ggtitle(x))
  dev.off()
}

celltype <- c("Exc-neuron",
              "Glial cell"
              # "Inh-neuron"
              )
for (i in celltype) {
  draw_scatter_merge_glial(i)
}
# calculate rescue ratio
non_sig_ratio_df <- matrix(ncol = 4)
non_sig_ratio_df <- as.data.frame(non_sig_ratio_df)
colnames(non_sig_ratio_df) <- c("name", "rescued_gene", "total_N1_N4_DEG", "ratio")
sig_ratio_df <- matrix(ncol = 4)
sig_ratio_df <- as.data.frame(sig_ratio_df)
colnames(sig_ratio_df) <- c("name", "rescued_DEG", "total_N1_N4_DEG", "ratio")

for (i in celltype) {
  # loading data
  N4_vs_N1 <- read.csv(paste0("./N4_vs_N1_re/N4_vs_N1_re_", i, ".csv"))
  overlap_df <- read.csv(paste0("./N4_N5_overlap_log2FC/N4_N5_DEG_overlap_", i, "_log2FC.csv"))
  # subset padj <= 0.05: only subset out SD3D vs. NS DEG
  N4_vs_N1_sig <- N4_vs_N1[N4_vs_N1$p_val_adj <= 0.05,]
  # extract rescued genes
  overlap_df_1 <- overlap_df[(overlap_df$avg_log2FC_N4_N1_DEGs>0 & overlap_df$avg_log2FC_N5_N4<0),]
  overlap_df_2 <- overlap_df[(overlap_df$avg_log2FC_N4_N1_DEGs<0 & overlap_df$avg_log2FC_N5_N4>0),]
  df <- rbind(overlap_df_1, overlap_df_2)
  non_sig_ratio_df[nrow(non_sig_ratio_df)+1,] <- c(i, nrow(df), nrow(N4_vs_N1_sig), nrow(df)/nrow(N4_vs_N1_sig))
  df <- df[df$N5_N4_sig == "sig",]
  sig_ratio_df[nrow(sig_ratio_df)+1,] <- c(i, nrow(df), nrow(N4_vs_N1_sig), nrow(df)/nrow(N4_vs_N1_sig))
}
write.csv(non_sig_ratio_df, file = "./N4_N5_overlap_log2FC/non_sig_ratio.csv")
write.csv(sig_ratio_df, file = "./N4_N5_overlap_log2FC/sig_ratio.csv")
