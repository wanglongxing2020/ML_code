library(tidyverse)
library(tidyr)
library(ggplot2)
library(cols4all)
library(ggpubr)
library(pheatmap)
library(ggvenn)
DAG <- read.csv("./DAG.marker.csv",header = T)
tcga <- read.csv("./LIHC_4529deg.csv",header = T,row.names = 1)
wgcna <- read.csv("./WGCNA_marker.csv",header = T)

DAG <- DAG[DAG$p_val_adj <0.05,]
#交集韦恩图
venn_list = list("scRNA-DAG-Marker"=DAG$X,"TCGA-Deg-Marker"=rownames(tcga),
                 "WGCNA-Marker"=wgcna$WGCNA_marker)
pdf('./venny_NEW.pdf',width = 10,height = 8)
ggvenn(venn_list,fill_color = c('#CD5C5C','#ADD8E6',"#DEB887"),fill_alpha = 0.5,
       set_name_size = 4,stroke_size = 0.7)
dev.off()


hubgene <-Reduce(intersect, list(DAG$X, rownames(tcga),
                                 wgcna$WGCNA_marker))


hubgene <- data.frame(hubgene)
write.csv(hubgene,file = "./37gene.csv")
