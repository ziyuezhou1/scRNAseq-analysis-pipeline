# split `so.renamed_res0.2.rds` by condition
# then run pyscenic for each condition
# loading packages
library(Seurat)

# loading data
seu <- readRDS("so.renamed_res0.2.rds")

# Split by sample
cells_by_group <- split(colnames(seu), seu$group_id)
seu <- lapply(cells_by_group, function(i) {
  subset(seu, cells = i)
})

for (i in 1:5) {
  saveRDS(seu[[i]], file = paste0("so.renamed_res0.2_N", i, ".rds"))
}

# save each group as loom
save_as_loom <- function(x){
  seu_x <- seu[[x]]
  exprMat <- seu_x@assays$RNA@data
  # write.table(exprMat, file = "./exprMat.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
  cellInfo <- seu_x@meta.data
  add_cell_annotation <- function(loom, cellAnnotation)
  {
    cellAnnotation <- data.frame(cellAnnotation)
    if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
    {
      warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
      cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
      cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
    }
    
    if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
    if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
    
    cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
    # Add annotation
    for(cn in colnames(cellAnnotation))
    {
      add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
    }
    
    invisible(loom)
  }
  
  loom <- build_loom(paste0("so.renamed_res.0.2_N", x, ".loom"), dgem=exprMat)
  loom <- add_cell_annotation(loom, cellInfo)
  close_loom(loom)
}

for (i in 1:5) {
  save_as_loom(i)
}
