# Differential composition analysis
# loading packages
library(scDC)
library(Seurat)
library(ggplot2)

# loading data
seu <- readRDS("../#3-Annotation-no-integration/3.2#-Annotation-res0.2/so.renamed_res0.2.rds")
png("./rename_clusters_dimplot_nolabel.png",
    width = 1200,
    height = 1200)
DimPlot(seu, reduction = "umap",label = FALSE, label.size = 10) + theme(aspect.ratio = 1)
dev.off()
sample <- seu$sample_id
celltypes <- Idents(seu)
group <- seu$group_id

table(sample)
t1 <- table(sample, celltypes)
write.csv(t1, "./celltype_number.csv")
info <- as.data.frame(table(group, celltypes))
info_N2_N1 <- subset(info, group %in% c('ZT6', 'ZT18'))
info_N3_N1 <- subset(info, group %in% c('ZT6', 'SD6h ZT0-ZT6'))
info_N3_N2 <- subset(info, group %in% c('ZT18', 'SD6h ZT0-ZT6'))
info_N4_N1 <- subset(info, group %in% c('ZT6', 'SD3D 8h/day (ZT0-ZT8)'))
info_N5_N1 <- subset(info, group %in% c('ZT6', 'SD3D R7'))

# compare N2 vs N1
info_N2_N1$group <- factor(info_N2_N1$group,
                           levels = c('ZT18', 'ZT6'))
png(filename = "celltype_proportion_N2_vs_N1.png",
    width = 800,
    height = 600)
ggplot(data = info_N2_N1,
       aes(x = celltypes, y = Freq, fill = group)) +
  geom_bar(position = 'fill', stat = 'identity') +
  coord_flip() +
  scale_fill_discrete(breaks = c('ZT6', 'ZT18'),
                      type = c('#03396c', '#b3cde0')) +
  theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0.50,
             linetype="dotted",
             color = "black",
             linewidth = 0.75)
dev.off()

# compare N3 vs N1
info_N3_N1$group <- factor(info_N3_N1$group,
                           levels = c('SD6h ZT0-ZT6', 'ZT6'))
png(filename = "celltype_proportion_N3_vs_N1.png",
    width = 1000,
    height = 600)
ggplot(data = info_N3_N1,
       aes(x = celltypes, y = Freq, fill = group)) +
  geom_bar(position = 'fill', stat = 'identity') +
  coord_flip() +
  scale_fill_discrete(breaks = c('ZT6', 'SD6h ZT0-ZT6'),
                      type = c('#03396c', '#b3cde0')) +
  theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0.50,
             linetype="dotted",
             color = "black",
             linewidth = 0.75)
dev.off()

# compare N3 vs N2
info_N3_N2$group <- factor(info_N3_N2$group,
                           levels = c('SD6h ZT0-ZT6', 'ZT18'))
png(filename = "celltype_proportion_N3_vs_N2.png",
    width = 1000,
    height = 600)
ggplot(data = info_N3_N2,
       aes(x = celltypes, y = Freq, fill = group)) +
  geom_bar(position = 'fill', stat = 'identity') +
  coord_flip() +
  scale_fill_discrete(breaks = c('ZT18', 'SD6h ZT0-ZT6'),
                      type = c('#03396c', '#b3cde0')) +
  theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0.50,
             linetype="dotted",
             color = "black",
             linewidth = 0.75)
dev.off()

# compare N4 vs N1
info_N4_N1$group <- factor(info_N4_N1$group,
                           levels = c('SD3D 8h/day (ZT0-ZT8)', 'ZT6'))
png(filename = "celltype_proportion_N4_vs_N1.png",
    width = 1000,
    height = 600)
ggplot(data = info_N4_N1,
       aes(x = celltypes, y = Freq, fill = group)) +
  geom_bar(position = 'fill', stat = 'identity') +
  coord_flip() +
  scale_fill_discrete(breaks = c('ZT6', 'SD3D 8h/day (ZT0-ZT8)'),
                      type = c('#03396c', '#b3cde0')) +
  theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0.50,
             linetype="dotted",
             color = "black",
             linewidth = 0.75)
dev.off()

# compare N5 vs N1
info_N5_N1$group <- factor(info_N5_N1$group,
                           levels = c('SD3D R7', 'ZT6'))
png(filename = "celltype_proportion_N5_vs_N1.png",
    width = 1000,
    height = 600)
ggplot(data = info_N5_N1,
       aes(x = celltypes, y = Freq, fill = group)) +
  geom_bar(position = 'fill', stat = 'identity') +
  coord_flip() +
  scale_fill_discrete(breaks = c('ZT6', 'SD3D R7'),
                      type = c('#03396c', '#b3cde0')) +
  theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0.50,
             linetype="dotted",
             color = "black",
             linewidth = 0.75)
dev.off()

# stacked bar chart
mycolor <- c("#47bb8a",
              "#ca66c8",
              "#92b540",
              "#6d71d8",
              "#cd9c2e",
              "#4f3788",
              "#5aa554",
              "#8c307d",
              "#33d4d1",
              "#d5474f",
              "#5f88d3",
              "#c4592f",
              "#b682d3",
              "#99933f",
              "#ca427f",
              "#c68845",
              "#db79b7",
              "#984629",
              "#892c5a",
              "#d96a66",
              "#94283d",
              "#de6b87",
              "#C9F2DC"
              )
png(filename = "./cell_ratio.png",
    width = 1000,
    height = 500)
ggplot(data = info,
       aes(x = group, y = Freq, fill = celltypes)) +
  geom_bar(position = 'fill', stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = mycolor) +
  theme(panel.background = element_blank(),
        text = element_text(size = 20))
dev.off()

# scDC analysis
res_scDC_noClust <- scDC_noClustering(celltypes, sample, calCI = TRUE,
                                      calCI_method = c("percentile", "BCa", "multinom"),
                                      nboot = 50)
png(filename = "./DCA.png",
    width = 3600,
    height = 600)
barplotCI(res_scDC_noClust, levels(group))
dev.off()
