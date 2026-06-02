#!/usr/bin/env Rscript
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
#library(ggunchull)
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
df <- df %>%
  mutate(log10p = -log10(p.adj),
         is_zero = (p.adj == 0))

max_finite <- max(df$log10p[is.finite(df$log10p)], na.rm = TRUE)
df$log10p[is.infinite(df$log10p)] <- max_finite * 1.2  
p1<-ggplot(df, aes(x = log10(FC), y = log10p)) +
  geom_point(aes(color = is_zero, shape = is_zero), size = 2) +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("p.adj > 0", "p.adj = 0")) +
  scale_shape_manual(values = c(16, 17), 
                     labels = c("p.adj > 0", "p.adj = 0")) +
  labs(x = "log10(Fold Change)", y = "-log10(adjusted p-value)",
       color = "Significance", shape = "Significance") +
  theme_bw() +
  geom_text_repel(data = subset(df, is_zero), aes(label = gene), 
                  nudge_y = 0.2, nudge_x = 0.1)  
p1 <- ggplot(df, aes(x = log10(FC), y = -log10(p.adj), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("#3B7EA1","#7A7A7A"), limits = c('sig', 'None')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93))+
  labs(x = 'beita', y = 'log10 p-value', color = '', title = '')+theme_bw() +theme(axis.line = element_line(colour = "black"))+theme(panel.border = element_blank())+ theme(panel.grid =element_blank())+ geom_hline(yintercept = 1.30103,linetype="dashed")

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
