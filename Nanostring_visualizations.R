renv::restore()

library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)

output_dir <- 'results'
dir.create(output_dir, recursive=TRUE)

### Heatmap for siADAR vs niNT Casp8 and MLKL knockout MEF experiment 
## DEG were determined using nSolver v4.0  default settings
# Raw and normalized data found in GSE200985

res1 <- read.csv("data/Data1.csv", header = T, sep = ",")
head(res1)
rownames(res1) <- res1[,1]
res1 <- subset(res1, select = -X)
head(res1)


### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

#threshold for siNT vs siADAR DKO
threshold1 <- res1$padj < padj.cutoff & abs(res1$log2FoldChange) > lfc.cutoff
length(which(threshold1))


res1$threshold1 <- threshold1
sig1 <- data.frame(subset(res1, threshold1==TRUE))


### Extract normalized expression for significant genes
sig1_ordered <- sig1[order(sig1$padj), ]
top20_sig1_genes <- rownames(sig1_ordered[1:20, ])
sig1_ordered_subset <- subset(sig1_ordered, select = -c(1,2,3,16,17,18,19,20))
norm_sig1 <- sig1_ordered_subset[rownames(sig1),]


### Set a color palette
heat.colors <- brewer.pal(11, "PiYG")

### Run pheatmap
pheatmap(norm_sig1, color = heat.colors, cluster_cols = F, cluster_rows = T, show_rownames=T,
         border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20,
         filename = paste0(output_dir,'/siADARvsniNT_heatmap.pdf'))

### Volcano plot for 'Casp8 KO ADAR siRNA vs Casp8 KO NT siRNA'
## DEG were determined using nSolver v4.0  default settings
# Raw and normalized data found in GSE200985

res2 <- read.csv("data/Data2.csv", header = T, sep = ",")
rownames(res2) <- res2[,1]
res2 <- subset(res2, select = -X)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

#threshold
threshold2 <- res2$P.value < padj.cutoff & abs(res2$Log2.fold.change) > lfc.cutoff
length(which(threshold2))
res2$threshold2 <- threshold2

# Volcano plot
res2_df <- data.frame(res2)

ggplot(res2_df) +
  geom_point(aes(x=Log2.fold.change, y=-log10(P.value), colour=threshold2)) +
  ggtitle("Casp8 KO ADAR siRNA vs Casp8 WT ADAR siRNA") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
ggsave(paste0(output_dir,'/Casp8KOADARsiRNAvsCasp8WTADARsiRNA.pdf'))

res2_df_ordered <- res2_df[order(res2_df$P.value), ] 

res2_df_ordered$genelabels <- rownames(res2_df_ordered) %in% rownames(res2_df_ordered[1:20,])

p1 <- ggplot(res2_df_ordered) +
  geom_point(aes(x = Log2.fold.change, y = -log10(P.value), colour = threshold2)) +
  geom_text_repel(aes(x = Log2.fold.change, y = -log10(P.value), label = ifelse(genelabels == T, rownames(res2_df_ordered),""))) +
  ggtitle("Casp8 KO ADAR siRNA vs Casp8 WT ADAR siRNA") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  xlim(-4.5,4.5) +
  ylim(0,6) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

p1 + theme_minimal()
ggsave(paste0(output_dir,'/Casp8KOADARsiRNAvsCasp8WTADARsiRNA_ordered.pdf'))

### WT vs ZBP1ko mutant ADAR1 MEF experiment
## DEG were determined using nSolver v4.0  default settings
# Raw and normalized data found in GSE200986

res3 <- read.csv("data/Data3.csv", header = T, sep = ",")
head(res3)
rownames(res3) <- res3[,1]
res3 <- subset(res3, select = -X)
head(res3)


### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

#threshold for siNT vs siADAR DKO 
threshold3 <- res3$padj <  lfc.cutoff
length(which(threshold3))

res3$threshold3 <- threshold3
sig3 <- data.frame(subset(res3, threshold3==TRUE))

### Set a color palette
heat.colors <- brewer.pal(11, "PiYG")

### Extract normalized expression for significant genes
sig3_ordered <- sig3[order(sig3$padj), ]
top20_sig3_genes <- rownames(sig3_ordered[1:20, ])
sig3_ordered_subset <- subset(sig3_ordered, select = -c(1,2,3,4,11,12,13,14,15))
norm_sig3 <- sig3_ordered_subset[rownames(sig3),]


### Run pheatmap
pheatmap(norm_sig3, color = heat.colors, cluster_cols = F, cluster_rows = T, show_rownames=T,
         border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20,
         filename = paste0(output_dir,'/WTvsZBP1ko.pdf'))
