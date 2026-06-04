#!/usr/bin/env Rscript
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidydr)
data = readRDS("E:\\MASLD_snRNA_seq_seurat_v4.rds")
meta<-scRNA@meta.data
scRNA$celltype.stim <- paste(scRNA$Cell_type_broad, scRNA$group, sep = "_")
seurat_obj<-subset(scRNA,celltype.stim==c("Hepatocyte_MASLD"))
seurat_obj
set.seed(42)
count_matrix <- GetAssayData(seurat_obj,assay="RNA",slot="counts")
target_gene_name="SUOX"
seurat_obj <- FindVariableFeatures(object=seurat_obj, selection.method="vst", nfeatures=10000)
high_variable_genes <- VariableFeatures(seurat_obj)
input_data <- as.data.frame(count_matrix[unique(c(target_gene_name, high_variable_genes)),])
library(scTenifoldNet.lite)
library(data.table)
ko_analysis_result <-scTenifoldKnk(countMatrix = input_data, gKO = "SUOX", qc_minLSize = 0, fast = TRUE)
diff_regulation_df <- ko_analysis_result$diffRegulation
#write.table(significant_diff_table, file="SUOX_KO.txt", sep="\t", quote=F, row.names=F)
df <- diff_regulation_df[diff_regulation_df$gene != c("SUOX"), ]
##KO
df[which(df$p.adj < 0.05),'sig'] <- 'Significant'
df[which(df$p.adj>= 0.05),'sig'] <- 'No Significant'
p1 <- ggplot(df, aes(x = log10(FC), y = -log10(p.adj), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("#3B7EA1","#7A7A7A"), limits = c('sig', 'None')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93))+
  labs(x = 'log10(FC)', y = '-log10(p.adj)', color = '', title = '')+theme_bw() +theme(axis.line = element_line(colour = "black"))+theme(panel.border = element_blank())+ theme(panel.grid =element_blank())+ geom_hline(yintercept = 1.30103,linetype="dashed")

p1
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
up <- df[df$gene%in%c("PLCG2", "NLRC5", "SERPINA1", "CALM1", "CFB","APOB","PRKCE","STK10"),]
p2 <- p1 + theme(legend.position = 'none') +
  geom_text_repel(data =up, aes(x = log10(FC), y = -log10(p.adj), label = gene),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = T)
p2

#Visualizing imputed gene expression
mydata<-read.table("magic_with_meta.tsv",header=T)
library(ggpointdensity)
ggplot(data=mydata, aes(x=mydata$SUOX, y=mydata$CTCF)) +
  geom_hex(bins = 80) +   
  scale_fill_viridis_c() 
cor.test(mydata$SUOX,mydata$CTCF,method="spearman")
