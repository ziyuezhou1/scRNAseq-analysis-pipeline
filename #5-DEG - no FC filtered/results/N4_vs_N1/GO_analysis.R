# loading package
library(DOSE)
library(org.Hs.eg.db)
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
library(biomaRt)
library(simplifyEnrichment)


# loading data
genes <- read.csv("./TFLink_targets_of_PBX3-homo.csv")
head(genes)

# GO
# convert ensembl into entrezid
entrez <- bitr(genes$Accession, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(entrez,2)
# BP
BP <- enrichGO(gene = entrez$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP", #choosing BP
                 pvalueCutoff = 0.05, #set p value cutoff
                 qvalueCutoff = 0.10,
                 readable = TRUE)
write.csv(summary(BP),file = "./TFLink_targets_of_PBX3_GOBP.csv",row.names =FALSE)
CC <- enrichGO(gene = entrez$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "CC", #choosing BP
               pvalueCutoff = 0.05, #set p value cutoff
               qvalueCutoff = 0.10,
               readable = TRUE)
write.csv(summary(CC),file = "./TFLink_targets_of_PBX3_GOCC.csv",row.names =FALSE)

# # BP_1 <- gofilter(BP, level = 3)
# # BP_1 <- BP_1[1:20,]
# print(BP)
mat = GO_similarity(BP$ID, ont = "BP")
BP_1 <- simplifyGO(mat)
# print(BP_1)
BP_1 <- merge(BP, BP_1, by.x = "ID", by.y = "id")
write.csv(BP_1,file = paste0("./",x,"_DEG_GOBP_clustered.csv"),row.names =FALSE)
  # BP_1 <- readxl::read_excel(paste0("./",x,"_DEG_GOBP_plot.xlsx"))
  # BP_1 <- BP_1[order(BP_1$p.adjust),]
  # png(filename = paste0("./",x,"_DEG_GOBP.png"),
  #     width = 1800,
  #     height = 1200)
  # p <- ggplot(BP_1,
  #        aes(y = fct_rev(reorder(Description, +p.adjust)), x = Count, fill = p.adjust))+
  #   geom_col(width = 0.8) +
  #   scale_x_continuous() +
  #   scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), trans = 'reverse')+
  #   labs(x = "Counts", y = NULL, title = paste0(x," DEG \nGO Biological Pathways"))+
  #   theme(text = element_text(size = 40), axis.text.y = element_text(size = 40))
  # print(p)
  # dev.off()


GO("nuc_up_p_0.06") # level4
