# loading packages


# circadian exclusive
celltype <- c("Astro",
              "Exc-neuron",
              "Inh-neuron",
              "Micro",
              "Oligo",
              "OPC")
for (i in celltype) {
  cir_exc <- read.csv(file = paste0("./circadian_genes/circadian_exclusive_DEGs_", i, ".csv"))
  cir_exc$avg_log2FC_N2_vs_N1 <- cir_exc$avg_log2FC.x
  cir_exc$avg_log2FC_N3_vs_N2 <- cir_exc$avg_log2FC.y
  cir_exc_N2_vs_N1_up <- cir_exc[cir_exc$avg_log2FC_N2_vs_N1 >= 0.585,]
  cir_exc_N2_vs_N1_down <- cir_exc[cir_exc$avg_log2FC_N2_vs_N1 <= -0.585,]
  cir_exc_N3_vs_N2_up <- cir_exc[cir_exc$avg_log2FC_N3_vs_N2 >= 0.585,]
  cir_exc_N3_vs_N2_down <- cir_exc[cir_exc$avg_log2FC_N3_vs_N2 <= -0.585,]
  cir_exc_res <- rbind(cir_exc_N2_vs_N1_up,
                       cir_exc_N2_vs_N1_down,
                       cir_exc_N3_vs_N2_up,
                       cir_exc_N3_vs_N2_down)
  write.csv(cir_exc_res, file = paste0("./circadian_genes/log2FC_filtered_circadian_exclusive_DEGs_", i, ".csv"))
}

# sleep exclusive DEGs
celltype <- c("Astro",
              "Exc-neuron",
              "Inh-neuron",
              "Micro",
              "Neuron",
              "Oligo")
for (i in celltype) {
  sleep_exc <- read.csv(file = paste0("./sleep_genes/sleep_exclusive_DEGs_", i, ".csv"))
  sleep_exc$avg_log2FC_N2_vs_N1 <- sleep_exc$avg_log2FC.x
  sleep_exc$avg_log2FC_N3_vs_N1 <- sleep_exc$avg_log2FC.y
  sleep_exc_N2_vs_N1_up <- sleep_exc[sleep_exc$avg_log2FC_N2_vs_N1 >= 0.585,]
  sleep_exc_N2_vs_N1_down <- sleep_exc[sleep_exc$avg_log2FC_N2_vs_N1 <= -0.585,]
  sleep_exc_N3_vs_N1_up <- sleep_exc[sleep_exc$avg_log2FC_N3_vs_N1 >= 0.585,]
  sleep_exc_N3_vs_N1_down <- sleep_exc[sleep_exc$avg_log2FC_N3_vs_N1 <= -0.585,]
  sleep_exc_res <- rbind(sleep_exc_N2_vs_N1_up,
                         sleep_exc_N2_vs_N1_down,
                         sleep_exc_N3_vs_N1_up,
                         sleep_exc_N3_vs_N1_down)
  write.csv(sleep_exc_res, file = paste0("./sleep_genes/log2FC_filtered_sleep_exclusive_DEGs_", i, ".csv"))
}

# R7 rescued DEGs
celltype <- c("Astro",
              "Exc-neuron",
              "Glial cell",
              "Inh-neuron",
              "Micro",
              "Oligo",
              "OPC")
for (i in celltype) {
  DEG <- read.csv(file = paste0("./N4_N5_overlap_log2FC/N4_N5_DEG_overlap_", i, "_log2FC.csv"))
  DEG_N5_N4_sig <- DEG[DEG$N5_N4_sig == "sig",]
  DEG_N4_N1_up <- DEG_N5_N4_sig[DEG_N5_N4_sig$avg_log2FC_N4_N1_DEGs >= 0.585,]
  DEG_N4_N1_down <- DEG_N5_N4_sig[DEG_N5_N4_sig$avg_log2FC_N4_N1_DEGs <= -0.585,]
  res <- rbind(DEG_N4_N1_up, DEG_N4_N1_down)
  write.csv(res, file = paste0("./acute_chronic_gene_analysis/log2FC_filtered_rescued_DEGs_", i, ".csv"))
}
