# loading packages
library(scuttle)
library(Seurat)
library(SingleCellExperiment)

# loading pseudobulk data
seurat <- readRDS("./data/so.renamed_res0.2.rds")
sce <- as.SingleCellExperiment(seurat)
agg.sce <- aggregateAcrossCells(sce, ids=sce$sample_id)
class(assay(agg.sce))
agg.sce <- calculateTPM(agg.sce)
agg.sce <- as.data.frame(agg.sce)
agg.sce$SYMBOL <- rownames(agg.sce)
# loading bulk data
bulk <- read.csv("./output_all-ExpDiff.csv")
genes <- bulk$Gene
genes <- sub('\\.[0-9]*$', '', genes)
bulk$Gene_without_dots <- genes
head(bulk)
symbol <- bitr(genes, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
bulk <- merge(bulk, symbol, by.x = "Gene_without_dots", by.y = "ENSEMBL")
head(bulk)
write.csv(bulk, "./nuc_bulk_output_all-ExpDiff_with_gene_names.csv")

# (2022). Machine Learning-Assisted Identification of Factors Contributing to the Technical Variability Between Bulk and Single-Cell RNA-Seq Experiments. e-pub ahead of print, doi: 10.21203/rs.3.rs-1247889/v1.
# overlap two dataframes
overlap <- merge(agg.sce, bulk, by = "SYMBOL")
nrow(overlap) ## 22030
overlap$NS_nuc_mean <- rowMeans(overlap[,c("NS_1_FC_clean_L001","NS_2_FC_clean_L002","NS_3_FC_clean_L003","NS_4_FC_clean_L004")])
overlap$SD1D_nuc_mean <- rowMeans(overlap[,c("SD1D_1_FC_clean_L005","SD1D_2_FC_clean_L006","SD1D_3_FC_clean_L007","SD1D_4_FC_clean_L008")])
overlap$SD3D_nuc_mean <- rowMeans(overlap[,c("SD3D_1_FC_clean_L009","SD3D_2_FC_clean_L010","SD3D_3_FC_clean_L011","SD3D_4_FC_clean_L012")])
overlap$R7_nuc_mean <- rowMeans(overlap[,c("R7_1_FC_clean_L013","R7_2_FC_clean_L014","R7_3_FC_clean_L015","R7_4_FC_clean_L016")])
overlap$NS_nuc_log2TPM <- log2(overlap$NS_nuc_mean + 0.001)
overlap$SD1D_nuc_log2TPM <- log2(overlap$SD1D_nuc_mean + 0.001)
overlap$SD3D_nuc_log2TPM <- log2(overlap$SD3D_nuc_mean + 0.001)
overlap$R7_nuc_log2TPM <- log2(overlap$R7_nuc_mean + 0.001)
overlap$NS_pseudo_log2TPM <- log2(overlap$N1 + 0.001)
overlap$SD1D_pseudo_log2TPM <- log2(overlap$N3 + 0.001)
overlap$SD3D_pseudo_log2TPM <- log2(overlap$N4 + 0.001)
overlap$R7_pseudo_log2TPM <- log2(overlap$N5 + 0.001)


# violin plot
data <- data.frame(
  group = c(rep("pseudo-bulk", nrow(overlap)), rep("nuc bulk", nrow(overlap))),
  log2TPM = c(overlap$NS_pseudo_log2TPM, overlap$NS_nuc_log2TPM)
)
head(data)
# Most basic violin chart
write.csv(summary(data[data$group == "nuc bulk",]), "./nuc_bulk_stat_summary_NS.csv")
write.csv(summary(data[data$group == "pseudo-bulk",]), "./nuc_pseudobulk_stat_summary_NS.csv")

png("./NS_nuc_bulk_vs_pseudobulk_violin.png")
p <- ggplot(data, aes(x=group, y=log2TPM, fill=group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() + geom_boxplot(width = 0.1) +
  labs(y = "log2(TPM + 0.001)", title = "NS nuc bulk vs. pseudobulk") +
  scale_fill_brewer(palette = "Blues") + theme_classic() +
  theme(text = element_text(size = 20))
print(p)
dev.off()

data <- data.frame(
  group = c(rep("pseudo-bulk", 4*nrow(overlap)), rep("nuc bulk", 4*nrow(overlap))),
  log2TPM = c(overlap$NS_pseudo_log2TPM,
              overlap$SD1D_pseudo_log2TPM,
              overlap$SD3D_pseudo_log2TPM,
              overlap$R7_pseudo_log2TPM,
              overlap$NS_nuc_log2TPM,
              overlap$SD1D_nuc_log2TPM,
              overlap$SD3D_nuc_log2TPM,
              overlap$R7_nuc_log2TPM)
)

write.csv(summary(data[data$group == "nuc bulk",]), "./bulk_stat_summary.csv")
write.csv(summary(data[data$group == "pseudo-bulk",]), "./pseudobulk_stat_summary.csv")

png("./nuc_bulk_vs_pseudobulk_violin.png")
p <- ggplot(data, aes(x=group, y=log2TPM, fill=group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() + geom_boxplot(width = 0.1) +
  labs(y = "log2(TPM + 0.001)", title = "nuc bulk vs. pseudo-bulk") +
  scale_fill_brewer(palette = "Blues") + theme_classic() +
  theme(text = element_text(size = 20))
print(p)
dev.off()

NS_NS_cor <- cor(overlap$NS_nuc_log2TPM, overlap$NS_pseudo_log2TPM, method = "pearson")
## 0.7929647
NS_SD1D_cor <- cor(overlap$NS_nuc_log2TPM, overlap$SD1D_pseudo_log2TPM, method = "pearson")
## 0.7982346
NS_SD3D_cor <- cor(overlap$NS_nuc_log2TPM, overlap$SD3D_pseudo_log2TPM, method = "pearson")
## 0.7950635
NS_R7_cor <- cor(overlap$NS_nuc_log2TPM, overlap$R7_pseudo_log2TPM, method = "pearson")
## 0.7912243
SD1D_SD1D_cor <- cor(overlap$SD1D_nuc_log2TPM, overlap$SD1D_pseudo_log2TPM, method = "pearson")
## 0.8014646
SD1D_SD3D_cor <- cor(overlap$SD1D_nuc_log2TPM, overlap$SD3D_pseudo_log2TPM, method = "pearson")
## 0.7993659
SD1D_R7_cor <- cor(overlap$SD1D_nuc_log2TPM, overlap$R7_pseudo_log2TPM, method = "pearson")
## 0.7950294
SD3D_SD3D_cor <- cor(overlap$SD3D_nuc_log2TPM, overlap$SD3D_pseudo_log2TPM, method = "pearson")
## 0.7986144
SD3D_R7_cor <- cor(overlap$SD3D_nuc_log2TPM, overlap$R7_pseudo_log2TPM, method = "pearson")
## 0.7943024
R7_R7_cor <- cor(overlap$R7_nuc_log2TPM, overlap$R7_pseudo_log2TPM, method = "pearson")
## 0.7940195