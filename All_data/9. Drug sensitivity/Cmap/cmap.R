degs = read.csv('./DNN_alldeg_mRNA(more).csv',header = T)  
up = degs[degs$logFC>0,]
down = degs[degs$logFC<0,]
up = up[order(up$logFC,decreasing = TRUE),]
up = up$X[1:150]

down = down[order(down$logFC,decreasing = FALSE),]
down = down$X[1:150]


write.csv(up,'up.csv')
write.csv(down,'down.csv')


#结果分析
result = read.delim('./LIHC_Cmap_export.txt')

res <- read.csv("./CmapScore_res.csv",header = T)
res$Drug <- factor(res$Drug,levels = res$Drug)
res$Database <- factor(res$Database,levels = c("CTRP","GDSC","PRSIM"))
library(ggplot2)
library(ggpubr)

pdf("./Cmap.pdf",width = 16,height = 8)
ggplot(data = res,mapping =aes(x=Drug,y=Cmap_score,color=Change,fill=Change))+geom_bar(size = 1.2,position="dodge", stat="identity",width = 0.8,alpha=0.7)+
  geom_text(aes(label=Cmap_score),position = position_dodge(0.9),color="black")+
  scale_color_manual(values=c("down"="#DD5F60","up"="#7DDFD7"))+
  scale_fill_manual(values=c("down"="#DD5F60","up"="#7DDFD7"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none")+
  facet_grid( ~Database,drop = T,scale="free",space = "free_x")+
  theme(axis.text=element_text(colour='black',size=11),
        strip.text.x = element_text(size = 20,colour = "white"),
        strip.background.x = element_rect(
          color="black", fill="#F2B379"),
        strip.placement = "outside",panel.spacing.x = unit(1, "cm"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

