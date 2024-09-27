library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(grid)
data <- read.csv("./ALL_hot_LIHC.csv",header = T,row.names = 1,check.names = F)

data2 <- data

for(i in 6:ncol(data2)){
  data2[,i] <-2*((data2[,i]-min(data2[,i]))/(max(data2[,i])-min(data2[,i])))-1
}

data3 <- t(data2)
data4 <- data3[-1,]
annotaion_col <- data2[,1,drop=F]
annotation_row <- read.csv("./Cell_type.csv",header = T,row.names = 1,check.names = F)
annotation_colors <- list(RS = colorRampPalette(c("navy","white","firebrick3"))(50),
                          Type= c(ESTIMATE="#E99C93",CIBERSORT="#9FACD3",
                                  MCPcounter="#F1DBB9",TIMER="#D9D1E3"))

coul <- colorRampPalette(brewer.pal(9, "OrRd"))(50)

pdf("./Cell_heatmap.pdf",width = 8,height = 16)
pheatmap(data4,annotation_col = annotaion_col,scale = "none",
         annotation_row = annotation_row,
         annotation_names_row = FALSE,
         annotation_colors = annotation_colors,
         show_colnames = F,cluster_rows = F,cluster_cols = F,
         color=coul,cellwidth=0.5, cellheight=10,
         gaps_row = c(4,26,36),border_color = NA)

dev.off()
