# loading package
library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(rlist)
library(VennDiagram)
library(tidyverse)
library(scales)

# loading data
N2_vs_N1 <- read.csv("../results/N2_vs_N1/N2_vs_N1_Inh-SST.csv")
# filtering
N2_vs_N1 <- N2_vs_N1[N2_vs_N1$p_val_adj <= 0.05,]
# SST GO
SST_DEG <- as.character(N2_vs_N1$X)
SSTGO <- bitr(SST_DEG, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
head(SSTGO,2)
# BP
SST_BP <- enrichGO(gene = SSTGO$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP", #choosing BP
                     pvalueCutoff = 0.05, #set p value cutoff
                     qvalueCutoff = 0.10,
                     readable = TRUE)
write.csv(summary(SST_BP),file = "./SST_GOBP.csv",row.names =FALSE)
SST_BP_1 <- SST_BP[order(SST_BP$p.adjust),]
SST_BP_1 <- SST_BP_1[1:20,]
png(filename = paste0("./SST_GOBP.png"),
    width = 1500,
    height = 1200)
ggplot(SST_BP_1,
       aes(y = fct_rev(reorder(Description, +p.adjust)), x = Count, fill = p.adjust))+
  geom_col(width = 0.8) +
  scale_x_continuous() +
  scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')+
  labs(x = "Counts", y = NULL, title = "SST NW vs. NS DEG \nGO Biological Pathways")+
  theme(text = element_text(size = 25), axis.text.y = element_text(size = 25))
# barplot(SST_BP,
#         showCategory=20,
#         title = paste0("SST NW vs. NS DEG \nGO Biological Pathways"),
#         font.size = 20)+
#   theme(plot.title = element_text(hjust = 0.5),
#         text = element_text(size = 20))+
#   scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
dev.off()

# upregulated DEG
# filtering
N2_vs_N1_up <- N2_vs_N1[N2_vs_N1$avg_log2FC > 0,]
# SST GO
SST_DEG_up <- as.character(N2_vs_N1_up$X)
SSTGO_up <- bitr(SST_DEG_up, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
head(SSTGO_up,2)
# BP
SST_BP_up <- enrichGO(gene = SSTGO_up$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP", #choosing BP
                   pvalueCutoff = 0.05, #set p value cutoff
                   qvalueCutoff = 0.10,
                   readable = TRUE)
write.csv(summary(SST_BP_up),file = "./SST_up_GOBP.csv",row.names =FALSE)
SST_BP_up_1 <- SST_BP_up[order(SST_BP_up$p.adjust),]
SST_BP_up_1 <- SST_BP_up_1[1:20,]
png(filename = paste0("./SST_up_GOBP.png"),
    width = 1500,
    height = 1200)
ggplot(SST_BP_up_1,
       aes(y = fct_rev(reorder(Description, +p.adjust)), x = Count, fill = p.adjust))+
  geom_col(width = 0.8) +
  scale_x_continuous() +
  scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')+
  labs(x = "Counts", y = NULL, title = "SST NW vs. NS upregulated DEG \nGO Biological Pathways")+
  theme(text = element_text(size = 25), axis.text.y = element_text(size = 25))
# barplot(SST_BP,
#         showCategory=20,
#         title = paste0("SST NW vs. NS DEG \nGO Biological Pathways"),
#         font.size = 20)+
#   theme(plot.title = element_text(hjust = 0.5),
#         text = element_text(size = 20))+
#   scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
dev.off()

# downregulated DEG
# filtering
N2_vs_N1_down <- N2_vs_N1[N2_vs_N1$avg_log2FC < 0,]
# SST GO
SST_DEG_down <- as.character(N2_vs_N1_down$X)
SSTGO_down <- bitr(SST_DEG_down, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
head(SSTGO_down,2)
# BP
SST_BP_down <- enrichGO(gene = SSTGO_down$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP", #choosing BP
                      pvalueCutoff = 0.05, #set p value cutoff
                      qvalueCutoff = 0.10,
                      readable = TRUE)
write.csv(summary(SST_BP_down),file = "./SST_down_GOBP.csv",row.names =FALSE)
SST_BP_down_1 <- SST_BP_down[order(SST_BP_down$p.adjust),]
SST_BP_down_1 <- SST_BP_down_1[1:20,]
png(filename = paste0("./SST_down_GOBP.png"),
    width = 1500,
    height = 1200)
ggplot(SST_BP_down_1,
       aes(y = fct_rev(reorder(Description, +p.adjust)), x = Count, fill = p.adjust))+
  geom_col(width = 0.8) +
  scale_x_continuous() +
  scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')+
  labs(x = "Counts", y = NULL, title = "SST NW vs. NS upregulated DEG \nGO Biological Pathways")+
  theme(text = element_text(size = 25), axis.text.y = element_text(size = 25))
# barplot(SST_BP,
#         showCategory=20,
#         title = paste0("SST NW vs. NS DEG \nGO Biological Pathways"),
#         font.size = 20)+
#   theme(plot.title = element_text(hjust = 0.5),
#         text = element_text(size = 20))+
#   scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')
dev.off()