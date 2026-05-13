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


GO_analysis_subtype <- function(x){
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
  write.csv(cirGO, file = paste0("./circadian_ENTREZID_", x, ".csv"))
  # BP
  cir_BP <- enrichGO(gene = cirGO$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP", #choosing BP
                     pvalueCutoff = 0.05, #set p value cutoff
                     qvalueCutoff = 0.10,
                     readable = TRUE)
  write.csv(summary(cir_BP),file = paste0("./circadian_genes/GO_circadian_DEG_", x, ".csv"),row.names =FALSE)
  BP_1 <- gofilter(cir_BP, level = 4)
  BP_1 <- BP_1[order(cir_BP$p.adjust),]
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

GO_analysis_subtype("Inh-SST")
