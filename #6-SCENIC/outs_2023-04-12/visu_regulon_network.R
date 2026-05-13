## loading packages
library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
library(SCENIC)
library(network)
library(igraph)

# define function
visuNetwork <- function(regulon.name){
  
  adj.ls <- regulon.name %>% map(~{
    tmp <- adj[which(adj$TF == .x), ]
    tmp <- tmp[order(tmp$importance, decreasing = T),]
    loci <- unique(c(1:50, grep('MUC2', tmp$target)))
    tmp <- tmp[loci,]
    return(tmp)
  })
  
  adj.sub <- Reduce('rbind', adj.ls)
  
  ## generate network
  edge.df <- adj.sub
  colnames(edge.df) <- c('from', 'to', 'weights')
  edge.df$from <- as.character(edge.df$from)
  edge.df$to <- as.character(edge.df$to)
  
  #net1 <- network(edge.df, directed = T, loops = T)
  #plot(net1)
  
  vertex <- unique(c(edge.df$from, edge.df$to))
  net <- graph_from_data_frame(d = edge.df, vertices = vertex, directed = T)
  
  
  ## vertex size scaled to log2FC
  #combined.tmp <- AverageExpression(combined, assays = 'RNA', features = vertex)
  #fc <- combined.tmp$RNA$SC_Hyper/combined.tmp$RNA$SC_AS
  #fc <- combined.tmp$RNA$SC_AS/combined.tmp$RNA$SC_Hyper
  #fc[which(fc > quantile(fc, 0.9))] <- quantile(fc, 0.9)
  #fc[which(fc < quantile(fc, 0.1))] <- quantile(fc, 0.1)
  #vsize <- scale(fc, center = min(fc), scale = max(fc) - min(fc))
  #vsize <- vsize+0.1
  
  ## vertex color: whether DEG in hyper vs as
  # vcl <- ifelse(vertex %in% dn.deg, 'blue', ifelse(vertex %in% up.deg, 'red', 'grey'))
  #  vlabel <- rep("", length(vertex))
  # vertex label 
  # loci <- which(vcl != 'grey' | vertex %in% edge.df$from | vertex == 'MUC2')
  #  vlabel[loci] = vertex[loci] 
  # vertex size scaled to weights
  vsize <- scale(edge.df$weights, center = min(quantile(edge.df$weights)), scale = max(quantile(edge.df$weights)) - min(quantile(edge.df$weights)))
  
  pdf(file = paste0(regulon.name, "_regulon_network.pdf"))
  plot(
    net,
    # vertex.label = vlabel,
    vertex.label.cex = 1,
    vertex.label.dist = 1,
    vertex.label.color = 'black',
    vertex.size = c(1, vsize+0.1)*10,
    #vertex.size = vsize*10,
    vertex.frame.color = NA,
    edge.arrow.size = 0.5,
    #  vertex.color = vcl, 
    main = regulon.name
  )
  dev.off()
  
  
}

# loading data
regulons_incidMat <- read.csv('so.renamed_res0.2_non_integrated_pySCENIC_regulons_incidence_matrix.csv')
regulons_incidMat_new <- regulons_incidMat[,-1]
rownames(regulons_incidMat_new) <- regulons_incidMat[,1]
regulons <- regulonsToGeneLists(regulons_incidMat_new)  # can convert to gene list for processing with R
head(regulons)
adj <- read.csv('pySCENIC_GRN_adjacencies.csv')
names(regulons) <- gsub("[_(+)]", "", names(regulons))
lapply(names(regulons)[grep('Erg', regulons)][1], visuNetwork)
