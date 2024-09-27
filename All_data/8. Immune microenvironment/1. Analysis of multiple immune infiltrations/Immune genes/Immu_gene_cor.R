library(tidyverse)
library(Hmisc)
immue <- read.csv("./immugene.csv",header = T)
load("~/TCGA_tumor_tpm.rda")

setdiff(immue$Gene,rownames(TCGA))
nogene <- c("IFNA1","IFNA2","IL2")
immue_gene <- immue$Gene[!immue$Gene%in%nogene]
TCGA <- TCGA[immue_gene,]

res<- t(TCGA) %>%as.data.frame()

#相关性分析
load("~/TCGA_RS.rda")
res2 <- res[rownames(TCGA_RS),]
dat <- TCGA_RS[,4,drop=F]

dat2 <- cbind(dat,res2)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res<-rcorr(as.matrix(dat2),type="pearson")

res1 <- flattenCorrMatrix(res$r, res$P)
res2 <- res1[res1$row=="RS",]
res2$p_value <- ifelse(res2$p < 0.05,ifelse(res2$p < 0.01,
                                            ifelse(res2$p < 0.001,
                                                   ifelse(res2$p < 0.0001,"****","***"),"**"),"*")
                       ,NA)
write.csv(res2,file = "Immue_cor_P.csv")

dat3 <- dat2[order(dat2$RS,decreasing = F),]

write.csv(dat3,file = "Immue_hot.csv")


#作图
data <- read.csv("./Hot_data.csv",header = T,row.names = 1,check.names = F)
data2 <- data

for(i in 2:ncol(data2)){
  data2[,i] <-2*((data2[,i]-min(data2[,i]))/(max(data2[,i])-min(data2[,i])))-1
}

data3 <- t(data2)
data4 <- data3[-1,]
annotaion_col <- data2[,1,drop=F]
annotation_row <- read.csv("./Hot_type.csv",header = T,row.names = 1,check.names = F)
annotation_colors <- list(RS = colorRampPalette(c("navy","white","firebrick3"))(50),
                          Type= c(`Antigen present`="#E99C93",`Cell adhesion`="#9FACD3",
                                  `Co-inhibitor`="#F1DBB9",`Co-stimulator`="#D9D1E3",
                                  Ligand ="#D0E2E8",other="#CCE2A2",Receptor="#C7AD8A"
                                  ))
coul <- colorRampPalette(brewer.pal(9, "OrRd"))(50)

pdf("./Immuegene_heatmap.pdf",width = 8,height = 10)
pheatmap(data4,annotation_col = annotaion_col,scale = "none",
         annotation_row = annotation_row,
         annotation_names_row = FALSE,
         annotation_colors = annotation_colors,
         show_colnames = F,cluster_rows = F,cluster_cols = F,
         color=coul,cellwidth=0.5, cellheight=10,
         gaps_row = c(13,16,23,26,45,51),border_color = NA)

dev.off()
