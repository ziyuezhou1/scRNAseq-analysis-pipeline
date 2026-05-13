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

# create result directory
if (!dir.exists("./acute_chronic_gene_analysis")) {dir.create("./acute_chronic_gene_analysis")}

acute_chronic <- function(x){
  N3_vs_N1 <- read.csv(file = paste0("./N3_vs_N1_re/N3_vs_N1_re_", x, ".csv"))
  N4_vs_N1 <- read.csv(file = paste0("./N4_vs_N1_re/N4_vs_N1_re_", x, ".csv"))
  # filtering
  N3_vs_N1 <- N3_vs_N1[N3_vs_N1$p_val_adj <= 0.05,]
  N4_vs_N1 <- N4_vs_N1[N4_vs_N1$p_val_adj <= 0.05,]
  # acute GO
  acute_DEG <- as.character(N3_vs_N1$X)
  acuteGO <- bitr(acute_DEG, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  head(acuteGO,2)
  # BP
  acute_BP <- enrichGO(gene = acuteGO$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP", #choosing BP
                     pvalueCutoff = 0.05, #set p value cutoff
                     qvalueCutoff = 0.10,
                     readable = TRUE)
  write.csv(summary(acute_BP),file = paste0("./acute_chronic_gene_analysis/GO_acute_DEG_", x, ".csv"),row.names =FALSE)
  png(filename = paste0("./acute_chronic_gene_analysis/GO_acute_DEG_", x, ".png"),
      width = 1200,
      height = 1200)
  p <- barplot(acute_BP, 
               showCategory=20,
               title = paste0("acute SD DEG in ", x, "\nGO Biological Pathways"))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
  print(p)
  dev.off()
  # chronic GO
  chronic_DEG <- as.character(N4_vs_N1$X)
  chronicGO <- bitr(chronic_DEG, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  head(chronicGO,2)
  # BP
  chronic_BP <- enrichGO(gene = chronicGO$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "BP", #choosing BP
                       pvalueCutoff = 0.05, #set p value cutoff
                       qvalueCutoff = 0.10,
                       readable = TRUE)
  write.csv(summary(chronic_BP),file = paste0("./acute_chronic_gene_analysis/GO_chronic_DEG_", x, ".csv"),row.names =FALSE)
  png(filename = paste0("./acute_chronic_gene_analysis/GO_chronic_DEG_", x, ".png"),
      width = 1200,
      height = 1200)
  p <- barplot(chronic_BP, 
               showCategory=20,
               title = paste0("chronic SD DEG in ", x, "\nGO Biological Pathways"))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
  print(p)
  dev.off()
  # draw venn graph to check overlap
  overlap <- list(N3_vs_N1 = N3_vs_N1$X,
                  N4_vs_N1 = N4_vs_N1$X)
  png(paste0("./acute_chronic_gene_analysis/acute_chronic_venn_", x, ".png"),
      width = 600,
      height = 600)
  p <- plot(euler(overlap), quantities = TRUE, main = paste0("acute_chronic_venn_", x))
  print(p)
  dev.off()
}

celltype <- c("Astro",
              "Endo",
              "Exc-neuron",
              "Glial cell(Egfr+)",
              "Inh-neuron",
              "Micro",
              "Oligo",
              "OPC")

for (i in celltype) {
  try(
    acute_chronic(i)
  )
}

# GO BP for all cell types
if (!dir.exists("./acute_chronic_gene_analysis/all")) {dir.create("./acute_chronic_gene_analysis/all")}
acute_chronic_all <- function(x){
  N3_vs_N1 <- read.csv(file = paste0("./N3_vs_N1/N3_vs_N1_", x, ".csv"))
  N4_vs_N1 <- read.csv(file = paste0("./N4_vs_N1/N4_vs_N1_", x, ".csv"))
  # filtering
  N3_vs_N1 <- N3_vs_N1[N3_vs_N1$p_val_adj <= 0.05,]
  N4_vs_N1 <- N4_vs_N1[N4_vs_N1$p_val_adj <= 0.05,]
  # acute GO
  acute_DEG <- as.character(N3_vs_N1$X)
  acuteGO <- bitr(acute_DEG, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  head(acuteGO,2)
  # BP
  acute_BP <- enrichGO(gene = acuteGO$ENTREZID,
                       OrgDb = org.Mm.eg.db,
                       ont = "BP", #choosing BP
                       pvalueCutoff = 0.05, #set p value cutoff
                       qvalueCutoff = 0.10,
                       readable = TRUE)
  write.csv(summary(acute_BP),file = paste0("./acute_chronic_gene_analysis/all/GO_acute_DEG_", x, ".csv"),row.names =FALSE)
  png(filename = paste0("./acute_chronic_gene_analysis/all/GO_acute_DEG_", x, ".png"),
      width = 1200,
      height = 1200)
  p <- barplot(acute_BP, 
               showCategory=20,
               title = paste0("acute SD DEG in ", x, "\nGO Biological Pathways"))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
  print(p)
  dev.off()
  # chronic GO
  chronic_DEG <- as.character(N4_vs_N1$X)
  chronicGO <- bitr(chronic_DEG, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  head(chronicGO,2)
  # BP
  chronic_BP <- enrichGO(gene = chronicGO$ENTREZID,
                         OrgDb = org.Mm.eg.db,
                         ont = "BP", #choosing BP
                         pvalueCutoff = 0.05, #set p value cutoff
                         qvalueCutoff = 0.10,
                         readable = TRUE)
  write.csv(summary(chronic_BP),file = paste0("./acute_chronic_gene_analysis/all/GO_chronic_DEG_", x, ".csv"),row.names =FALSE)
  png(filename = paste0("./acute_chronic_gene_analysis/all/GO_chronic_DEG_", x, ".png"),
      width = 1200,
      height = 1200)
  p <- barplot(chronic_BP, 
               showCategory=20,
               title = paste0("chronic SD DEG in ", x, "\nGO Biological Pathways"))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
  print(p)
  dev.off()
  # draw venn graph to check overlap
  overlap <- list(N3_vs_N1 = N3_vs_N1$X,
                  N4_vs_N1 = N4_vs_N1$X)
  png(paste0("./acute_chronic_gene_analysis/all/acute_chronic_venn_", x, ".png"),
      width = 600,
      height = 600)
  p <- plot(euler(overlap), quantities = TRUE, main = paste0("acute_chronic_venn_", x))
  print(p)
  dev.off()
}

acute_chronic_all("Inh-PV&SST")

# compare GO in different major cell types
all <- data.frame()
for (i in celltype) {
  DEG1 <- read.csv(file = paste0("./N3_vs_N1_re/N3_vs_N1_re_", i, ".csv"))
  DEG2 <- read.csv(file = paste0("./N4_vs_N1_re/N4_vs_N1_re_", i ,".csv"))
  DEG1 <- DEG1[DEG1$p_val_adj <= 0.05,]
  DEG2 <- DEG2[DEG2$p_val_adj <= 0.05,]
  DEG1$SYMBOL <- DEG1$X
  DEG2$SYMBOL <- DEG2$X
  find1 <- bitr(DEG1$X, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  find2 <- bitr(DEG2$X, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  DEG1 <- merge(DEG1, find1, by = "SYMBOL")
  DEG2 <- merge(DEG2, find2, by = "SYMBOL")
  DEG1$group <- i
  DEG2$group <- i
  DEG1$othergroup <- "acute"
  DEG2$othergroup <- "chronic"
  all_DEG <- rbind(DEG1, DEG2)
  all <- rbind(all, all_DEG)
}
# GO: BP
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                fun = "enrichGO",
                                ont = "BP",
                                data = all,
                                OrgDb = "org.Mm.eg.db"
                                )
re <- gofilter(re, level = 3)
png("./acute_chronic_gene_analysis/all_compare_GOBP.png",
    width = 1500,
    height = 1200)
dotplot(re, x = "group") +
  facet_grid(~othergroup) +
  theme(text = element_text(size = 20))
dev.off()
write.csv(summary(re), file = "./acute_chronic_gene_analysis/acute_chronic_all_BP_new_20230222.csv")
# GO: MF
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                      fun = "enrichGO",
                                      ont = "MF",
                                      data = all,
                                      OrgDb = "org.Mm.eg.db"
)
re <- gofilter(re, level = 3)
png("./acute_chronic_gene_analysis/all_compare_GOMF.png",
    width = 1500,
    height = 1200)
dotplot(re, x = "group") +
  facet_grid(~othergroup) +
  theme(text = element_text(size = 20))
dev.off()
write.csv(summary(re), file = "./acute_chronic_gene_analysis/acute_chronic_all_MF_new_20230222.csv")
# GO: CC
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                      fun = "enrichGO",
                                      ont = "CC",
                                      data = all,
                                      OrgDb = "org.Mm.eg.db"
)
re <- gofilter(re, level = 3)
png("./acute_chronic_gene_analysis/all_compare_GOCC.png",
    width = 1500,
    height = 1200)
dotplot(re, x = "group") +
  facet_grid(~othergroup) +
  theme(text = element_text(size = 20))
dev.off()
write.csv(summary(re), file = "./acute_chronic_gene_analysis/acute_chronic_all_CC_new_20230222.csv")
# KEGG
re <- clusterProfiler::compareCluster(geneClusters = ENTREZID~group+othergroup,
                                      fun = "enrichKEGG",
                                      data = all,
                                      organism = "mmu"
)
png("./acute_chronic_gene_analysis/all_compare_KEGG.png",
    width = 1500,
    height = 1200)
dotplot(re, x = "group") +
  facet_grid(~othergroup) +
  theme(text = element_text(size = 20))
dev.off()
write.csv(summary(re), file = "./acute_chronic_gene_analysis/acute_chronic_all_KEGG_new_20230222.csv")
