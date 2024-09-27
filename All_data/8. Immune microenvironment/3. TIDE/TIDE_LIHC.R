library(ggpubr)
library(ggplot2)
library(tidyverse)
library(ggsci)
load("~/TCGA_RS.rda")
rt <- TCGA_RS[order(TCGA_RS$RS, decreasing = F),]
rt$RS_group <- ifelse(rt$RS_group=="high","High","Low")
rt$id <- c(1:length(rt$RS)) # 新增id列  
rt$id2 <- paste(rt$RS_group, rt$id, sep = '_') # 将风险分组和id串联 
colnames(rt)[6] <- "Patient"



TIDE <- read.csv("./TIDE_LIHC_output.csv",header = T)
TIDE <- TIDE[,c(1,4)]

T_C <- merge(TIDE,rt,by="Patient")

T_C$RS_group <- factor(T_C$RS_group,levels = c("High","Low"))

pdf("./TIDE_deg.pdf",width = 6,height = 6)
ggplot(T_C,aes(x=RS_group,y=TIDE))+
  geom_point(position = 'jitter',aes(color=RS_group),size=2,alpha=1)+theme_bw()+
  geom_violin(aes(fill=RS_group),color='grey',scale = 'width',linewidth=0.8,trim=T,alpha=0.7)+
  scale_color_jama()+
  scale_fill_jama()+
  stat_compare_means(comparisons = list(c("High","Low")),
                     method = "wilcox.test")+
  geom_boxplot(aes(color=RS_group),width=0.2,size=0.8,alpha=0.7)+
  xlab("MPCDI")+ylab("TIDE score")+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank()) 
dev.off()


pdf("./TIDE_cor.pdf",width = 6,height = 6)
ggplot(T_C,aes(x=RS,y=TIDE))+
  geom_point(color="#800000")+
  geom_smooth(color="#FFA500",fill="#FFA500",method = "lm",alpha=0.2)+
  stat_cor(method = "pearson")+theme_bw()+
  xlab("MPCDI")+ylab("TIDE score")
dev.off()
