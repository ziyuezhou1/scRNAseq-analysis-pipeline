# loading package
library(graphics)

# creating barplot matrix
# list celltype
celltype <- c("Astro",
              "Endo",
              "Exc-neuron",
              "Glial cell(Egfr+)",
              "Inh-neuron",
              "Micro",
              "Oligo",
              "OPC")
group <- c("N3 vs. N1", "N4 vs. N1")
bar_mat <- matrix(nrow = 2, ncol = 0)
rownames(bar_mat) <- group

for (i in celltype) {
  # loading data
  N3_vs_N1_i <- read.csv(file = paste0("./N3_vs_N1_re/N3_vs_N1_re_",i,".csv"))
  N4_vs_N1_i <- read.csv(file = paste0("./N4_vs_N1_re/N4_vs_N1_re_",i,".csv"))
  # filter out only sig
  N3_vs_N1_i_sig <- N3_vs_N1_i[N3_vs_N1_i$p_val_adj <= 0.05,]
  N4_vs_N1_i_sig <- N4_vs_N1_i[N4_vs_N1_i$p_val_adj <= 0.05,]
  tmp <- matrix(nrow = 2, ncol = 1)
  rownames(tmp) <- group
  colnames(tmp) <- i
  tmp[1,1] <- nrow(N3_vs_N1_i_sig)
  tmp[2,1] <- nrow(N4_vs_N1_i_sig)
  # combine into result matrix
  bar_mat <- cbind(bar_mat, tmp)
}

# draw bar plot
png("./barplot_N4_N5_DEG.png",
    width = 1300,
    height = 1500)
barplot(bar_mat,
        main = "DEG number in N3 vs. N1 and N4 vs. N1",
        xlab = "DEG number",
        ylab = "Cell Type",
        col = c("#2a8cd1", "#16537e"),
        beside = TRUE,
        horiz = TRUE,
        legend = TRUE,
        border = NA,
        cex.main = 2,
        cex.names = 2,
        cex.axis = 2,
        cex.lab = 1.5)
legend("topright",
       legend = rownames(bar_mat),
       pch = 15,
       col = c("#2a8cd1", "#16537e"),
       cex = 2)
dev.off()