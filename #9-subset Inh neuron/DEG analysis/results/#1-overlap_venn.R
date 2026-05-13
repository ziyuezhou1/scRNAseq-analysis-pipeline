# loading package
library(BioVenn)
library(eulerr)
library(UpSetR)
library(ggplot2)
# creating directory
if (!dir.exists("./N1_N2_N3_venn/")) {dir.create("./N1_N2_N3_venn/")}
if (!dir.exists("./N1_N4_N5_DEG_venn/")) {dir.create("./N1_N4_N5_DEG_venn/")}

draw_venn_3 <- function(x){
  # loading data
  N2_vs_N1 <- read.csv(paste0("./N2_vs_N1/N2_vs_N1_",x,".csv"))
  N3_vs_N1 <- read.csv(paste0("./N3_vs_N1/N3_vs_N1_",x,".csv"))
  N3_vs_N2 <- read.csv(paste0("./N3_vs_N2/N3_vs_N2_",x,".csv"))
  
  # subset significant rows
  N2_vs_N1 <- N2_vs_N1[N2_vs_N1$p_val_adj <= 0.05,]
  N3_vs_N1 <- N3_vs_N1[N3_vs_N1$p_val_adj <= 0.05,]
  N3_vs_N2 <- N3_vs_N2[N3_vs_N2$p_val_adj <= 0.05,]
  N2_vs_N1_genes <- as.list(N2_vs_N1$X)
  N3_vs_N1_genes <- as.list(N3_vs_N1$X)
  N3_vs_N2_genes <- as.list(N3_vs_N2$X)
  # draw venn graph
  mylist <- list(NW_vs_NS = N2_vs_N1_genes,
                 SD1D_vs_NS = N3_vs_N1_genes,
                 SD1D_vs_NW = N3_vs_N2_genes)
  png(paste0("./N1_N2_N3_venn/",x, "_N1_N2_N3.png"),
      width = 600,
      height = 600)
  p <- plot(euler(mylist), quantities = TRUE, main = paste0(x, "_NS_NW_SD1D"), cex = 1.5,
            cex.lab = c(10, 10, 10, 10, 10, 10))
  print(p)
  dev.off()
  
  png(paste0("./N1_N2_N3_venn/",x, "_N1_N2_N3_1.png"),
      width = 600,
      height = 600)
  p <- upset(fromList(mylist), nsets = 3,
           point.size = 5, line.size = 2,
           text.scale = c(3,3,1.3,1.3,2,3))
  print(p)
  dev.off()
}

draw_venn_2 <- function(x){
  # loading data
  N4_vs_N1 <- read.csv(paste0("./N4_vs_N1/N4_vs_N1_",x,".csv"))
  N5_vs_N4 <- read.csv(paste0("./N5_vs_N4/N5_vs_N4_",x,".csv"))

  # subset significant rows
  N4_vs_N1 <- N4_vs_N1[N4_vs_N1$p_val_adj <= 0.05,]
  N5_vs_N4 <- N5_vs_N4[N5_vs_N4$p_val_adj <= 0.05,]
  N4_vs_N1_genes <- as.vector(N4_vs_N1$X)
  N5_vs_N4_genes <- as.vector(N5_vs_N4$X)
  mylist <- list(N4_vs_N1 = N4_vs_N1_genes,
                 N5_vs_N4 = N5_vs_N4_genes)
  # draw venn graph
  png(paste0("./N1_N4_N5_DEG_venn/",x, "_N1_N4_N5_DEG.png"),
      width = 600,
      height = 600)
  p <- plot(euler(mylist), quantities = TRUE, main = paste0(x, "_N1_N4_N5_DEG"),
            cex.main = 3,
            cex.sub = 3,
            cex.lab = 3,
            cex.axis = 3)
  print(p)
  dev.off()
}

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


for (i in celltype) {
  # draw N1, N2, N3 DEG
  try(draw_venn_3(i))
  # draw N4, N5 DEG
  try(draw_venn_2(i))
}

