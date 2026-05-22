##Code for the scRNA-seq analysis
##KO
df<-read.table("E:\\2026\\MASLD\\SUOX_KO.txt",header=T)
df[which(df$p.adj < 0.05),'sig'] <- 'sig'
df[which(df$p.adj>= 0.05),'sig'] <- 'None'

# 计算 -log10(p.adj)，处理无穷大
# 处理 -log10(p.adj)
df <- df %>%
  mutate(log10p = -log10(p.adj),
         is_zero = (p.adj == 0))

max_finite <- max(df$log10p[is.finite(df$log10p)], na.rm = TRUE)
df$log10p[is.infinite(df$log10p)] <- max_finite * 1.2  # 将 p.adj=0 的点放在顶部

# 检查行数
print(nrow(df))  # 应为 10

# 绘图：X 轴使用 log10(FC)
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
                  nudge_y = 0.2, nudge_x = 0.1)  # 仅标记 p.adj=0 的基因
##
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

#
p2 <- p1 + theme(legend.position = 'none') +
  geom_text_repel(data =up, aes(x = log10(FC), y = -log10(p.adj), label = gene),
                  size = 5,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = T)
p2
c("#2C7BB6", "#D7191C", "#999999")
p2<-p1
p2+p1

##Vislizat in R
mydata<-read.table("magic_with_meta_raw2.tsv",header=T)
library(ggpointdensity)
ggplot(data=mydata, aes(x=mydata$PPP5C, y=mydata$ETS1)) +
  geom_hex(bins = 80) +   # bins 控制六边形数量
  scale_fill_viridis_c()  # 颜色方案


