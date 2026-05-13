# loading package
library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(rlist)
library(VennDiagram)
library(eulerr)
library(tidyverse)
library(ggvenn)
library(UpSetR)

# create result directory
if (!dir.exists("./circadian_genes")) {dir.create("./circadian_genes")}
if (!dir.exists("./sleep_genes")) {dir.create("./sleep_genes")}

GO_analysis <- function(x){
  # loading data
  N2_vs_N1 <- read.csv(paste0("./N2_vs_N1/N2_vs_N1_", x, ".csv"))
  N3_vs_N1 <- read.csv(paste0("./N3_vs_N1/N3_vs_N1_", x, ".csv"))
  N3_vs_N2 <- read.csv(paste0("./N3_vs_N2/N3_vs_N2_", x, ".csv"))
  # subset padj <= 0.05
  N2_vs_N1 <- N2_vs_N1[N2_vs_N1$p_val_adj <= 0.05,]
  N3_vs_N1 <- N3_vs_N1[N3_vs_N1$p_val_adj <= 0.05,]
  N3_vs_N2 <- N3_vs_N2[N3_vs_N2$p_val_adj <= 0.05,]
  # find overlap genes: possible rescued DEGs
  overlap_cir <- as.vector(intersect(N2_vs_N1$X, N3_vs_N2$X))
  overlap_sleep <- as.vector(intersect(N2_vs_N1$X, N3_vs_N1$X))
  # subset other data out
  N2_vs_N1_cir <- N2_vs_N1[N2_vs_N1$X %in% overlap_cir,]
  N2_vs_N1_sleep <- N2_vs_N1[N2_vs_N1$X %in% overlap_sleep,]
  N3_vs_N2_cir <- N3_vs_N2[N3_vs_N2$X %in% overlap_cir,]
  N3_vs_N1_sleep <- N3_vs_N1[N3_vs_N1$X %in% overlap_sleep,]
  # merge
  circadian <- merge(N2_vs_N1_cir,
                     N3_vs_N2_cir,
                     by = "X")
  sleep <- merge(N2_vs_N1_sleep,
                 N3_vs_N1_sleep,
                 by = "X")
  # export genes
  write.csv(circadian, file = paste0("./circadian_genes/circadian_DEGs_", x, ".csv"))
  write.csv(sleep, file = paste0("./sleep_genes/sleep_DEGs_", x, ".csv"))
  # GO analysis
  keytypes(org.Mm.eg.db)
  # circadian DEGs
  cirGO <- bitr(overlap_cir, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  head(cirGO,2)
  # BP
  cir_BP <- enrichGO(gene = cirGO$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP", #choosing BP
                     pvalueCutoff = 0.05, #set p value cutoff
                     qvalueCutoff = 0.10,
                     readable = TRUE)
  write.csv(summary(cir_BP),file = paste0("./circadian_genes/GO_circadian_DEG_", x, ".csv"),row.names =FALSE)
  BP_1 <- cir_BP[order(cir_BP$p.adjust),]
  BP_1 <- BP_1[1:20,]
  png(filename = paste0("./circadian_genes/GO_circadian_DEG_", x, ".png"),
      width = 1200,
      height = 1200)
  # p <- barplot(cir_BP, 
  #         showCategory=20,
  #         title = paste0("circadian DEG in ", x, "\nGO Biological Pathways"))+
  #   theme(plot.title = element_text(hjust = 0.5))+
  #   scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
  p <- ggplot(BP_1,
         aes(y = fct_rev(reorder(Description, +p.adjust)), x = Count, fill = p.adjust))+
    geom_col(width = 0.8) +
    scale_x_continuous() +
    scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')+
    labs(x = "Counts", y = NULL, title = paste0("circadian DEG in ", x, "\nGO Biological Pathways"))+
    theme(text = element_text(size = 25), axis.text.y = element_text(size = 25))
  print(p)
  dev.off()
  # sleep DEGs
  sleepGO <- bitr(overlap_sleep, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  head(sleepGO,2)
  # BP
  sleep_BP <- enrichGO(gene = sleepGO$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "BP", #choosing BP
                       pvalueCutoff = 0.05, #set p value cutoff
                       qvalueCutoff = 0.10,
                       readable = TRUE)
  write.csv(summary(sleep_BP),file = paste0("./sleep_genes/GO_sleep_DEG_", x, ".csv"),row.names =FALSE)
  png(filename = paste0("./sleep_genes/GO_sleep_DEG_", x, ".png"),
      width = 1200,
      height = 1200)
  p <- barplot(sleep_BP, 
          showCategory=20,
          title = paste0("sleep DEG in ", x, "\nGO Biological Pathways"))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
  print(p)
  dev.off()
}

celltype <- c("Inh-PV",
              "Inh-SST(Calb1+)",
              "Inh-VIP(Penk+)",
              "Inh-Lamp5(Dbp+)",
              "Inh-Adarb2+Lhx6-Pbx3+",
              "Inh-Erbb4+",
              "Inh-PAX6",
              "Inh-SST(Npy+Adgrg6+)")
for (i in celltype) {
  try(
    GO_analysis(i)
  )
}

# draw overlap venn
# circadian
celltype <- c(
  "Inh-PV",
  "Inh-SST(Calb1+)",
  "Inh-VIP(Penk+)",
  "Inh-Lamp5(Dbp+)",
  "Inh-Adarb2+Lhx6-Pbx3+",
  "Inh-Erbb4+",
  "Inh-PAX6",
  "Inh-SST(Npy+Adgrg6+)"
)
cir_list <- list()
for (i in celltype) {
  cir_DEG <- read.csv(file = paste0("./circadian_genes/circadian_DEGs_", i, ".csv"))
  cir_DEG <- cir_DEG$X
  cir_list[[i]] <- cir_DEG
}
png("./circadian_genes/circadian_upset.png",
    width = 1500,
    height = 800)
plot.new()
p <- upset(fromList(cir_list), nsets = 9,
           point.size = 5, line.size = 2,
           text.scale = c(3,3,1.3,1.3,2,3))
print(p)
dev.off()
# png(paste0("./circadian_genes/circadian_venn_all_celltypes_1.png"),
#     width = 600,
#     height = 600)
# p <- plot(euler(cir_list), quantities = TRUE, main = "circadian_venn_all_celltypes")
# # p <- ggvenn(cir_list)
# print(p)
# dev.off()

# sleep
celltype <- c(
  "Inh-PV",
  "Inh-SST(Calb1+)",
  "Inh-VIP(Penk+)",
  "Inh-Lamp5(Dbp+)",
  "Inh-Adarb2+Lhx6-Pbx3+",
  "Inh-Erbb4+",
  "Inh-PAX6",
  "Inh-SST(Npy+Adgrg6+)"
)
sleep_list <- list()
for (i in celltype) {
  sleep_DEG <- read.csv(file = paste0("./sleep_genes/sleep_DEGs_", i, ".csv"))
  sleep_DEG <- sleep_DEG$X
  sleep_list[[i]] <- sleep_DEG
}
png("./sleep_genes/sleep_upset.png",
    width = 1500,
    height = 800)
p <- upset(fromList(sleep_list), nsets = 9,
           point.size = 5, line.size = 2,
           text.scale = c(3,3,1.3,1.3,2,3))
print(p)
dev.off()
# png(paste0("./sleep_genes/sleep_venn_all_celltypes.png"),
#     width = 600,
#     height = 600)
# p <- plot(euler(sleep_list), quantities = TRUE, main = "sleep_venn_all_celltypes")
# print(p)
# dev.off()

# compare GO in different cell types
all <- data.frame()
for (i in celltype) {
  cir_DEG <- read.csv(paste0("./circadian_genes/circadian_DEGs_", i, ".csv"))
  sleep_DEG <- read.csv(paste0("./sleep_genes/sleep_DEGs_", i, ".csv"))
  cir_DEG$SYMBOL <- cir_DEG$X
  sleep_DEG$SYMBOL <- sleep_DEG$X
  cir_DEG <- cir_DEG[,c("SYMBOL", "avg_log2FC.x", "avg_log2FC.y", "p_val_adj.x", "p_val_adj.y")]
  sleep_DEG <- sleep_DEG[,c("SYMBOL", "avg_log2FC.x", "avg_log2FC.y", "p_val_adj.x", "p_val_adj.y")]
  if (nrow(cir_DEG) > 0 & nrow(sleep_DEG > 0)) {
    find_cir <- bitr(cir_DEG$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
    find_sleep <- bitr(sleep_DEG$SYMBOL, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
    cir_DEG <- merge(cir_DEG, find_cir, by = "SYMBOL")
    sleep_DEG <- merge(sleep_DEG, find_sleep, by = "SYMBOL")
    cir_DEG$group <- i
    sleep_DEG$group <- i
    cir_DEG$othergroup <- "circadian"
    sleep_DEG$othergroup <- "sleep"
    all_DEG <- rbind(cir_DEG, sleep_DEG)
    all <- rbind(all, all_DEG)} else {
      if (nrow(sleep_DEG) == 0 || nrow(cir_DEG) == 0) {next}
    }
}
all$group <- factor(all$group, levels = unique(all$group))
all$group <- factor(all$group, levels = c("Inh-PV",
                                          "Inh-SST(Calb1+)",
                                          "Inh-VIP(Penk+)",
                                          "Inh-Lamp5(Dbp+)",
                                          "Inh-Adarb2+Lhx6-Pbx3+",
                                          "Inh-Erbb4+",
                                          "Inh-PAX6",
                                          "Inh-SST(Npy+Adgrg6+)"))

# GO: BP
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                      fun = "enrichGO",
                                      ont = "BP",
                                      data = all,
                                      OrgDb = "org.Mm.eg.db")

re <- gofilter(re, level = 4)
re@compareClusterResult$group <- factor(re@compareClusterResult$group, levels = unique(re@compareClusterResult$group))
png("./circadian_sleep_BP_all_compare.png",
    width = 2600,
    height = 1400)
ggplot(re, 
       aes(y = fct_rev(reorder(Description, +p.adjust)),
            x = group)) +
  facet_grid(~othergroup) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_size_continuous(range = c(1, 20)) +
  scale_color_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse') +
  theme_bw() +
  labs(x = "Cell Type", y = NULL, title = NULL)+
  theme(text = element_text(size = 28),
        axis.text.y = element_text(size = 28),
        axis.text.x = element_text(size = 25),
        strip.text.x = element_text(size = 30))
# dotplot(re, x = "group") +
#   facet_grid(~othergroup) +
#   theme(text = element_text(size = 25),
#         axis.text.x.bottom = element_text(size = 25),
#         axis.text.y.left = element_text(size = 25))
dev.off()
write.csv(summary(re), file = "./circadian_sleep_BP_all_compare.csv")
# GO: MF
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                      fun = "enrichGO",
                                      ont = "MF",
                                      data = all,
                                      OrgDb = "org.Mm.eg.db")
re <- gofilter(re, level = 4)
re@compareClusterResult$group <- factor(re@compareClusterResult$group, levels = unique(re@compareClusterResult$group))
png("./circadian_sleep_MF_all_compare.png",
    width = 1500,
    height = 1200)
dotplot(re, x = "group") +
  facet_grid(~othergroup) +
  theme(text = element_text(size = 20))
dev.off()
write.csv(summary(re), file = "./circadian_sleep_MF_all_compare.csv")
# GO: CC
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                      fun = "enrichGO",
                                      ont = "CC",
                                      data = all,
                                      OrgDb = "org.Mm.eg.db")
re <- gofilter(re, level = 4)
re@compareClusterResult$group <- factor(re@compareClusterResult$group, levels = unique(re@compareClusterResult$group))
png("./circadian_sleep_CC_all_compare.png",
    width = 1500,
    height = 1200)
dotplot(re, x = "group") +
  facet_grid(~othergroup) +
  theme(text = element_text(size = 20))
dev.off()
write.csv(summary(re), file = "./circadian_sleep_CC_all_compare.csv")
# KEGG
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                      fun = "enrichKEGG",
                                      data = all,
                                      organism = "mmu")
re@compareClusterResult$group <- factor(re@compareClusterResult$group, levels = unique(re@compareClusterResult$group))
png("./circadian_sleep_KEGG_all_compare.png",
    width = 1500,
    height = 1200)
dotplot(re, x = "group") +
  facet_grid(~othergroup) +
  theme(text = element_text(size = 20))
dev.off()
write.csv(summary(re), file = "./circadian_sleep_KEGG_all_compare.csv")

# check core genes in circadian genes
cir_celltype <- c("Inh-PV",
                  "Inh-SST(Calb1+)",
                  "Inh-VIP(Penk+)",
                  "Inh-Lamp5(Dbp+)",
                  "Inh-Adarb2+Lhx6-Pbx3+",
                  "Inh-Erbb4+",
                  "Inh-PAX6",
                  "Inh-SST(Npy+Adgrg6+)")
cir_list <- list()
for (i in cir_celltype) {
  cir <- read.csv(file = paste0("./circadian_genes/circadian_DEGs_", i, ".csv"))
  cir_list[[i]] <- cir
}
cir_join <- Reduce(function(x, y) merge(x, y, by = "X"), cir_list)
write.csv(cir_join, file = "./circadian_genes/circadian_core_genes.csv")

# find exclusive circadian genes in each celltype
cir_celltype <- c("Inh-PV",
                  "Inh-SST(Calb1+)",
                  "Inh-VIP(Penk+)",
                  "Inh-Lamp5(Dbp+)",
                  "Inh-Adarb2+Lhx6-Pbx3+",
                  "Inh-Erbb4+",
                  "Inh-PAX6",
                  "Inh-SST(Npy+Adgrg6+)")
all_cir_list <- list()
for (i in cir_celltype) {
  all_cir_list[[i]] <- read.csv(file = paste0("./circadian_genes/circadian_DEGs_", i, ".csv"))
}
find_exclusive_genes <- function(x){
  want <- all_cir_list[[x]]
  rest_cir_list <- all_cir_list
  rest_cir_list[[x]] <- NULL
  for (i in rest_cir_list) {
    want <- subset(want, !(X %in% i$X))
  }
  write.csv(want, file = paste0("./circadian_genes/circadian_exclusive_DEGs_", x, ".csv"), row.names = FALSE)
}
for (j in cir_celltype) {
  find_exclusive_genes(j)
}

# find exclusive sleep genes in each celltype
sleep_celltype <- c(
  "Inh-PV",
  "Inh-SST(Calb1+)",
  "Inh-VIP(Penk+)",
  "Inh-Lamp5(Dbp+)",
  "Inh-Adarb2+Lhx6-Pbx3+",
  "Inh-Erbb4+",
  "Inh-PAX6",
  "Inh-SST(Npy+Adgrg6+)")
all_sleep_list <- list()
for (i in sleep_celltype) {
  all_sleep_list[[i]] <- read.csv(file = paste0("./sleep_genes/sleep_DEGs_", i, ".csv"))
}
find_exclusive_genes <- function(x){
  want <- all_sleep_list[[x]]
  rest_sleep_list <- all_sleep_list
  rest_sleep_list[[x]] <- NULL
  for (i in rest_sleep_list) {
    want <- subset(want, !(X %in% i$X))
  }
  write.csv(want, file = paste0("./sleep_genes/sleep_exclusive_DEGs_", x, ".csv"), row.names = FALSE)
}
for (j in sleep_celltype) {
  find_exclusive_genes(j)
}

