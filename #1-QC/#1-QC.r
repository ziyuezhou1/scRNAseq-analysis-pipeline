#QC
library(Seurat)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scds)

# LOAD AND REFORMAT DATA
# Load raw counts
fastq_dirs <- list.dirs("./data/", recursive = FALSE, full.names = TRUE)
names(fastq_dirs) <- basename(fastq_dirs)
sce <- DropletUtils::read10xCounts(fastq_dirs)
# 文件名要对应地改成barcodes，genes，matrix，否则无法导入
sce

# Rename row/colData colnames & sce dimnames
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL", "Gene Expression")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))
dimnames(sce) <- list(
  with(rowData(sce), paste(SYMBOL, sep = ".")),
  with(colData(sce), paste(sample_id, barcode, sep = "."))
)

# Load metadata
md_dir <- file.path("F:\\dataset20230202/", "dataset_metadata.xlsx")
md <- readxl::read_excel(md_dir)
m <- match(sce$sample_id, md$`Sample Name`)
sce$group_id <- md$sorting_condition[m]

dim(sce)

# DOUBLET REMOVAL
# Split by sample
cs_by_s <- split(colnames(sce), sce$sample_id)
sce_by_s <- lapply(cs_by_s, function(cs) sce[, cs])

# Run 'scds'
sce_by_s <- lapply(sce_by_s, function(u) {
  cxds_bcds_hybrid(bcds(cxds(u)))
})

# Remove doublets
sce_by_s <- lapply(sce_by_s, function(u) {
  # Compute expected nb. of doublets (10x)
  n_dbl <- ceiling(0.01 * ncol(u)^2 / 1e3)
  # Remove 'n_dbl' cells w/ highest doublet score
  o <- order(u$hybrid_score, decreasing = TRUE)
  u[, -o[seq_len(n_dbl)]]
})

# Merge back into single SCE
sce <- do.call(cbind, sce_by_s)
dim(sce)

# Create seurat object
seu <- CreateSeuratObject(counts = counts(sce),
                          project = "SDdataset",
                          meta.data = data.frame(colData(sce)),
                          min.cells = 3,
                          min.features = 200,)
seu



# use PercentageFeatureSet function to calculate mitochondrion QC
seu[['percent.mt']] <- PercentageFeatureSet(seu, pattern = "^mt-")
# Violin plot visualization
png("./vlnplot_qc_nol3.png",
    width = 1000,
    height = 800)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# FeatureScatter: visualization of relationship between two patterns
plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("./Featurescatter.png",
    width = 1200,
    height = 600)
plot1 + plot2
dev.off()
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8599307/
# In general, for mouse tissues, the 5% threshold performs well to distinguish between healthy and low-quality cells. 
# cutoff for nFeature_RNA is subjective
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8750 & percent.mt < 5)
dim(seu)
saveRDS(seu, file = "./seu.rds")