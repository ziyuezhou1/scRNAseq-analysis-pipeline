# circadian DEGs shared in neurons
# loading data
Exc <- read.csv(file = "./circadian_genes/circadian_DEGs_Exc-neuron.csv")
Inh <- read.csv(file = "./circadian_genes/circadian_DEGs_Inh-neuron.csv")

Neu <- merge(x = Exc,
             y = Inh,
             by = "X")
write.csv(Neu, file = "./circadian_genes/circadian_DEGs_Neuron.csv")

# circadian DEGs shared in glial cells
cir_celltype <- c("Astro",
                  "Micro",
                  "Oligo")
cir_list <- list()
for (i in cir_celltype) {
  cir <- read.csv(file = paste0("./circadian_genes/circadian_DEGs_", i, ".csv"))
  cir_list[[i]] <- cir
}
cir_join <- Reduce(function(x, y) merge(x, y, by = "X"), cir_list)
write.csv(cir_join, file = "./circadian_genes/circadian_DEGs_in_glia.csv")
