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
library(biomaRt)
library(simplifyEnrichment)
library(stringr)
library(scales)

GO <- function(x){
  BP_1 <- readxl::read_excel(paste0("./",x,"_DEG_GOBP_plot.xlsx"))
  # BP_1$`Function Cluster` <- factor(BP_1$`Function Cluster`, levels = unique(BP_1$`Function Cluster`))
  # BP_1$Description <- str_wrap(BP_1$Description, width = 60)
  png(filename = paste0("./",x,"_DEG_GOBP.png"),
      width = 2500,
      height = 1800)
  p <- ggplot(BP_1,
              aes(y = fct_rev(reorder(Description, p.adjust)), x = Count, fill = p.adjust))+
    geom_col(width = 0.8) +
    # facet_grid(`Function Cluster` ~ .,
    #            scales = "free",
    #            space = "free",
    #            labeller = labeller(`Function Cluster` = label_wrap_gen(width = 5))) +
    scale_x_continuous() +
    scale_fill_gradientn(colours = c("#005582","#9f45b0","#cc3300"), space = "Lab")+
    labs(x = "Counts", y = NULL, title = paste0(x," circadian DEG \nGO Biological Pathways"))+
    scale_y_discrete(labels = wrap_format(40)) +
    theme(text = element_text(size = 50),
          axis.text.y = element_text(size = 65),
          legend.key.height = unit(2.5,"cm"))
  print(p)
  dev.off()
}
GO("Inh-SST")