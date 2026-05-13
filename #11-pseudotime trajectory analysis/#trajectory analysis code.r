#trajectory analysis code
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(scran)
library(Seurat)
library(SingleCellExperiment)
library(dyno)
library(tidyverse)
so <- readRDS(file.path("../#5-DEG - no FC filtered/data/", "so.renamed_res0.2.rds"))
#split by sample
so <- subset(x = so, idents = c("Inh-PV&SST","Exc-L2&L4/5IT","Exc-neuron1","Exc-L6CT",
                                "Exc-neuron2","Inh-VIP&Lamp5","Exc-neuron3","Exc-L5/6N",
                                "Exc-neuron4","Interneuron","Exc-L2/3IT","Exc-Car3"), invert = FALSE)
# only Inh neuron
# so <- subset(x = so, idents = c("Inh-PV&SST","Inh-VIP&Lamp5","Interneuron"), invert = FALSE)

cells_by_sample <- split(colnames(so), so$sample_id)
so <- lapply(cells_by_sample, function(i) {
  subset(so, cells = i)
})
#analyze N1
sce <- as.SingleCellExperiment(so$N1, assay = "RNA")
dataset <- wrap_expression(
  expression = t(logcounts(sce)),
  counts = t(counts(sce))
)
dataset <- add_grouping(
  dataset, so$N1@active.ident
)
#select methods
guidelines <- guidelines_shiny(dataset)
# Reproduces the guidelines as created in the shiny app
# for exc and inh neurons
# Reproduces the guidelines as created in the shiny app
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE, 
#   expect_topology = FALSE, 
#   expected_topology = NULL, 
#   expect_cycles = TRUE, 
#   n_cells = 7505, 
#   n_features = 34867, 
#   memory = "6GB", 
#   prior_information = c("start_id", "end_id", "end_n", "start_n", "leaves_n", "groups_n", "features_id", "dimred"), 
#   docker = TRUE
# )
# guidelines <- dynguidelines::guidelines(answers = answers)
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE,
#   expect_topology = FALSE,
#   expected_topology = NULL,
#   expect_cycles = TRUE,
#   n_cells = 1498,
#   n_features = 34867,
#   memory = "10GB",
#   prior_information = c("start_id", "end_id", "end_n", "start_n", "leaves_n", "groups_n", "features_id", "dimred"),
#   docker = TRUE
# )
# guidelines <- dynguidelines::guidelines(answers = answers)
#come back to Rstudio after selecting
dataset <- add_prior_information(dataset, start_id = "N1.GGCTTGGTCCAGCCTT-1")
methods_selected <- guidelines$methods_selected
methods_selected
#Running the methods
library(akima)
library(hdf5r)
model <- infer_trajectory(dataset, dplyr::first(methods_selected))
pdf("./N1_trajectory.pdf",
    width = 10,
    height = 8)
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
dev.off()
pdf("./N1_trajectory_pseudotime.pdf",
    width = 10,
    height = 8)
plot_dimred(
  model,
  "pseudotime",
  pseudotime = calculate_pseudotime(model))+ 
  ggtitle("Pseudotime")
dev.off()
pdf("./N1_trajectory_heatmap.pdf",
    width = 10,
    height = 8)
plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  features_oi = 50
)
dev.off()

#analyze N4
sce <- as.SingleCellExperiment(so$N4, assay = "RNA")
dataset <- wrap_expression(
  expression = t(logcounts(sce)),
  counts = t(counts(sce))
)
dataset <- add_grouping(
  dataset, so$N4@active.ident
)
#select methods
guidelines <- guidelines_shiny(dataset)
# Reproduces the guidelines as created in the shiny app
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE, 
#   expect_topology = FALSE, 
#   expected_topology = NULL, 
#   expect_cycles = TRUE, 
#   n_cells = 6700, 
#   n_features = 34867, 
#   memory = "6GB", 
#   prior_information = c("start_id", "end_id", "end_n", "start_n", "leaves_n", "groups_n", "features_id", "dimred"), 
#   docker = TRUE
# )
# guidelines <- dynguidelines::guidelines(answers = answers)
# Reproduces the guidelines as created in the shiny app
# answers <- dynguidelines::answer_questions(
#   multiple_disconnected = FALSE,
#   expect_topology = FALSE,
#   expected_topology = NULL,
#   expect_cycles = TRUE,
#   n_cells = 2112,
#   n_features = 34867,
#   memory = "10GB",
#   prior_information = c("start_id", "end_id", "end_n", "start_n", "leaves_n", "groups_n", "features_id", "dimred"),
#   docker = TRUE
# )
# guidelines <- dynguidelines::guidelines(answers = answers)
#come back to Rstudio after selecting
dataset <- add_prior_information(dataset, start_id = "N4.CTTTCGGGTTATAGAG-1")
methods_selected <- guidelines$methods_selected
methods_selected
#Running the methods
library(akima)
model <- infer_trajectory(dataset, dplyr::first(methods_selected))
library(hdf5r)
pdf("./N4_trajectory.pdf",
    width = 10,
    height = 8)
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)
dev.off()
#pseudotime
pdf("./N4_pseudotime.pdf",
    width = 10,
    height = 8)
plot_dimred(
  model,
  "pseudotime",
  pseudotime = calculate_pseudotime(model))+
  ggtitle("Pseudotime")
dev.off()
#topfeature
pdf("./N4_heatmap.pdf",
    width = 10,
    height = 8)
plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  features_oi = 50
)
dev.off()
#branch_point_gene
branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()
branch_feature_importance <- calculate_branching_point_feature_importance(model, expression_source=dataset$expression, milestones_oi = branching_milestone)
branching_point_features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)
plot_heatmap(
  model,
  expression_source = dataset$expression,
  features_oi = branching_point_features
)
space <- dyndimred::dimred_mds(dataset$expression)

map(branching_point_features[1:12], function(feature_oi) {
  plot_dimred(model, dimred = space, expression_source = dataset$expression, feature_oi = feature_oi, label_milestones = FALSE) +
    theme(legend.position = "none") +
    ggtitle(feature_oi)
}) %>% patchwork::wrap_plots()               