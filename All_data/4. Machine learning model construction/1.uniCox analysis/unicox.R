library(survival)
library(survminer)
library(tidyverse)
library(glmnet)

load("../TCGA_tumor_tpm.rda")
clin <- read.csv("../clin_all.csv",header = T,row.names = 1)
hubgene <- read.csv("../37gene.csv",header = T,row.names = 1)
hubgene <- hubgene$hubgene
data <- TCGA[rownames(TCGA)%in%hubgene,]

data <- t(data) %>%as.data.frame()
clin1 <- clin[,c(1,2,3)]
rownames(clin1) <- clin1$Sample
clin1 <- clin1[,-1]
clin1$time <- clin1$time/365

rt <- merge(clin1,data,by="row.names")
rownames(rt) <- rt$Row.names
rt <- rt[,-1]
rt[,3:ncol(rt)] <-apply(rt[,3:ncol(rt)],2,as.numeric)

save(rt,file = "./TCGA_exp_clin39.rda")



#单因素cox分析
outTab=data.frame()
for (i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(time,status) ~ rt[,i],data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  coxP[is.na(coxP)] <- 1
  #保留预后相关的基因
  if(coxP<0.05){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.csv(outTab,file = "./TCGA_outTab_16.csv",row.names = F)


#-----绘图
data <- outTab
data$id <- factor(data$id,levels = rev(data$id))
data[2:5] <- apply(data[2:5],2,as.numeric)
data$PN <- ifelse(data$HR > 1 ,"P","N")
data$co <- ifelse(data$PN=="P","#ce5f06","#95c9e5")
p <- ggplot(data,aes(color = PN,x=HR,xmin=HR.95L,xmax = HR.95H,y=id))+
  #geom_pointrange(size = 0,fatten = 0)+   ##代码中有两处geom_pointrange，此处的代码参数设置为0，是为了避免报错和使背景色块位于误差棒下面
  scale_color_manual(values = c("#95c9e5","#ce5f06"),labels =c("Protective (P<0.05)","Risky (P<0.05)"))+
  theme_bw()+
  scale_x_continuous(limits = c(0.5,2.1),breaks = c(0.5,1,1.5,2))+
  xlab("HR(95%CI)")+ #X轴标题
  ylab("")+
  theme(
    legend.title = element_blank(),
    axis.text.y = element_text(color = rev(data$co)), #设置Y轴字体颜色
    axis.ticks.x = element_blank(), # 删除X轴刻度
    axis.ticks.y = element_blank() # 删除Y轴刻度
  )+
  geom_pointrange(
    size=4, # 控制线的宽度
    fatten = 1)+  #点的大小
  geom_vline(xintercept = 1, color = "black",size = 1,linetype="dotted")+ #添加直线
  guides(color=guide_legend(override.aes = list(size=1.5),reverse = T))

pdf("forest_16.pdf",width = 10,height = 12)
p
dev.off()


