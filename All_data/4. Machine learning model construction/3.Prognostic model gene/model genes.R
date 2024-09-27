library(ggplot2)
library(ggpubr)
library(tidyverse) 
library(ggforce)
library(gghalves)
library(ggdist) 

load("~/LIHC_expr_tpm.rda")

gene <- read.csv("./gene_4.csv",header = T)
data <- LIHC_expr
data <- data[gene$gene,]
data <- t(data) %>%as.data.frame()
data <- log2(data+1)
rownames(data) <- str_sub(rownames(data),1,16)
group <- sapply(strsplit(rownames(data),"\\-"),'[',4)
group <- sapply(strsplit(group,""),'[',1)
group <- gsub("2","1",group)
group <- ifelse(as.numeric(group) == 0,"Tumor","Normal")
data<-cbind(group,data)
data[,2:ncol(data)] <-apply(data[,2:ncol(data)] ,2,as.numeric)



Normaldat<- data[data$group=="Normal",]
Tumordat<- data[data$group=="Tumor",]

#新加
load("~/TCGA_tumor_tpm.rda")
Tumordat <- Tumordat[substr(rownames(Tumordat),1,12)%in%colnames(TCGA),]
Tumordat <- Tumordat[substr(rownames(Tumordat),15,15)==1,]
data <- rbind(Normaldat,Tumordat)

for(i in 2:ncol(data)){
  dat <- data[,c(1,i)]
  gene <- colnames(data)[i]
  boxplot <- ggplot(dat,aes(x=group,y=get(gene),fill=group))+
      geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
      geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
      geom_jitter(aes(fill=group),shape=21,size=2.5,width=0.2)+
      scale_fill_manual(values = c("#5cc3e8","#ffdb00"))+
      theme_bw()+
      stat_compare_means(comparisons = list(c("Normal","Tumor")),
                     method = "wilcox.test",
                     label = "p.signif")+
      theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=15),
        legend.position = "none")+
      ylab(paste0(gene," expression"))+xlab("")
  #输出图片
  pdf(file = paste0(gene, ".pdf"),width = 5.5,height = 5)
  print(boxplot)
  dev.off()
}


