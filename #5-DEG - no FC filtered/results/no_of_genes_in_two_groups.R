# function to count gene numbers
genenumber <- function(x, y){
  data <- read.csv(paste0("./", y, "_re/", y, "_re_", x, ".csv"))
  gene <- data[!(data$pct.1 == 0 & data$pct.2 ==0),]
  print(y)
  print(x)
  print(nrow(gene))
}

celltype <- c(
  "Astro",
  "Endo",
  "Exc-neuron",
  "Glial cell",
  "Inh-neuron",
  "Micro",
  "Neuron",
  "Oligo",
  "OPC"
)
for (i in celltype) {
  genenumber(i, "N2_vs_N1")
}
for (i in celltype) {
  genenumber(i, "N3_vs_N1")
}
for (i in celltype) {
  genenumber(i, "N4_vs_N1")
}
for (i in celltype) {
  genenumber(i, "N5_vs_N1")
}
