library(tidyverse)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(tidyr)
#TCGA-LIHC
load("~/TCGA_tumor_tpm.rda")
clin <- read.csv("./clin_all.csv",header = T,row.names = 1)
data <- t(TCGA) %>%as.data.frame()
clin1 <- clin[,c(1,2,3)]
rownames(clin1) <- clin1$Sample
clin1 <- clin1[,-1]
clin1$time <- clin1$time/365

rt <- merge(clin1,data,by="row.names")
rownames(rt) <- rt$Row.names
rt <- rt[,-1]
rt[,3:ncol(rt)] <-apply(rt[,3:ncol(rt)],2,as.numeric)

save(rt,file = "./TCGA_exp_clinALL.rda")

ALL_Cindex <- read.csv("./All_Cindex.csv",header = T,row.names = 1)
#MPCDI
TCGA_MPCDI_cindex <- ALL_Cindex[1,]
TCGA_MPCDI_cindex$cohort <- "MPCDI"
#1.Xiao-Wei Fu  
Coef <- c(0.0182,0.0005,0.0188)
#0.0182*GSDME+0.0005*GPX4+0.0188*SCAF11
gene <- c("GSDME","GPX4","SCAF11")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"GSDME"]*Coef[1]+dat1[,"GPX4"]*Coef[2]+dat1[,"SCAF11"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6154
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_1_cindex <- data.frame("cohort"='Xiao-Wei Fu',"C-index"=cindex["C"],
                              "Lower"=lower_cindex,"Upper"=upper_cindex)

#2.Jongmin Kim 
Coef <- c(-0.333,-0.4,0.339,0.387)
gene <- c("CDH1","ID2","MMP9","TCF3")
intersect(gene,colnames(rt)) #都在
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CDH1"]*Coef[1]+dat1[,"ID2"]*Coef[2]+
  dat1[,"MMP9"]*Coef[3]+dat1[,"TCF3"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6489
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_2_cindex <- data.frame("cohort"='Jongmin Kim',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

#3.Zhen Zhang
Coef <- c(-0.13681,-0.07452,-0.09026,-0.1212)
gene <- c("CAT","EHHADH","ALDH5A1","SLC27A5")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CAT"]*Coef[1]+dat1[,"EHHADH"]*Coef[2]+dat1[,"ALDH5A1"]*Coef[3]+
  dat1[,"SLC27A5"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6082
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_3_cindex <- data.frame("cohort"='Zhen Zhang',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

#4.Yue-ling Peng
Coef <- c(0.0029,0.0086,0.0024,0.9751)
gene <- c("HSP90AA1","PPIA","SQSTM1","USP21")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"HSP90AA1"]*Coef[1]+dat1[,"PPIA"]*Coef[2]+dat1[,"SQSTM1"]*Coef[3]+
  dat1[,"USP21"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6061
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_4_cindex <- data.frame("cohort"='Yue-ling Peng',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

#5.Xiangkun Wang
Coef <- c(-0.5,-0.592,-0.616)
gene <- c("JAK2","STAT5A","STAT6")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"JAK2"]*Coef[1]+dat1[,"STAT5A"]*Coef[2]+
  dat1[,"STAT6"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.657
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_5_cindex <- data.frame("cohort"='Xiangkun Wang',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#6.Shanshan Lu
Coef <- c(0.267389,0.489842)
gene <- c("ANO10", "CLCN2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ANO10"]*Coef[1]+dat1[,"CLCN2"]*Coef[2]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.64
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_6_cindex <- data.frame("cohort"='Shanshan Lu',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#7.Jia-heng Xing
Coef <- c(0.0310975867942454,0.0385177274200257,0.014946834324991)
gene <- c("AKT1", "MAPK3" ,"CASP3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"AKT1"]*Coef[1]+dat1[,"MAPK3"]*Coef[2]+
  dat1[,"CASP3"]*Coef[3]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6557
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_7_cindex <- data.frame("cohort"='Jia-heng Xing',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#8.Hong Peng
Coef<-c(0.0634,0.12544,0.0481,0.0483)
gene <- c("ATP7A", "MTF1", "GLS", "CDKN2A")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ATP7A"]*Coef[1]+dat1[,"MTF1"]*Coef[2]+
  dat1[,"GLS"]*Coef[3]+dat1[,"CDKN2A"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6535
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_8_cindex <- data.frame("cohort"='Hong Peng',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#9.Huaxiang Wang
Coef<-c(0.0618,0.0538,0.1552,0.0428,-0.3199,0.1208)
gene <- c("BUB3", "IGF2BP3","RBM3", "ILF3", "ZC3H13", "CCT3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"BUB3"]*Coef[1]+dat1[,"IGF2BP3"]*Coef[2]+
  dat1[,"RBM3"]*Coef[3]+dat1[,"ILF3"]*Coef[4]+
  dat1[,"ZC3H13"]*Coef[5]+dat1[,"CCT3"]*Coef[6]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6467
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_9_cindex <- data.frame("cohort"='Huaxiang Wang',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#10.Yuqin Tang
Coef<-c(1.3,0.68,1.17,1.49,1.17)
gene <- c("PSRC1", "SOCS2", "TMEM45A", "CCT5", "STC2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"PSRC1"]*Coef[1]+dat1[,"SOCS2"]*Coef[2]+
  dat1[,"TMEM45A"]*Coef[3]+dat1[,"CCT5"]*Coef[4]+
  dat1[,"STC2"]*Coef[5]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6425
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_10_cindex <- data.frame("cohort"='Yuqin Tang',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

TCGA_all_Cindex <- bind_rows(list(TCGA_MPCDI_cindex,TCGA_1_cindex,TCGA_2_cindex
                                   ,TCGA_3_cindex,TCGA_4_cindex,TCGA_5_cindex,
                                   TCGA_6_cindex,TCGA_7_cindex,TCGA_8_cindex,
                                   TCGA_9_cindex,TCGA_10_cindex))
TCGA_all_Cindex <-TCGA_all_Cindex[order(TCGA_all_Cindex$C.index,decreasing = T),]
write.csv(TCGA_all_Cindex,file = "./TCGA_all_Cindex.csv",row.names = F)

#Plot---------------
data <- read.csv("./TCGA_all_Cindex.csv",header = T)
data$cohort <- factor(data$cohort,levels = rev(data$cohort))
dat1 <- gather(data,C,value,1:3)
pdf(file = "./TCGA-LIHC.pdf",width = 5,height = 4)
ggplot(dat1,aes(cohort,value)) +
  #scale_y_continuous(breaks  = c(0.5,0.6,0.7),limits = c(0.5,0.75))+
  geom_segment(data=data,aes(x=cohort,y=Lower,xend=cohort,yend=Upper),
               color="#175676",size=1,alpha=0.5)+
  geom_point(data=subset(dat1,C=="C.index"),aes(color=C),size=6)+
  scale_color_manual(values = c("firebrick"))+
  geom_hline(aes(yintercept=0.6),linetype=5,col="#175676",alpha=0.5)+
  coord_flip()+
  labs(title = "TCGA-LIHC",x="",y="")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none")
dev.off()

#GSE14520
load("~/GSE14520_exp_clin.rda")
rt <- GSE14520_exp_clin[,-c(1:3)]
ALL_Cindex <- read.csv("./All_Cindex.csv",header = T,row.names = 1)
#MPCDI
GSE14520_MPCDI_cindex <- ALL_Cindex[2,]
GSE14520_MPCDI_cindex$cohort <- "MPCDI"
#1.Xiao-Wei Fu  
Coef <- c(0.0182,0.0005,0.0188)
#0.0182*GSDME+0.0005*GPX4+0.0188*SCAF11
gene <- c("DFNA5","GPX4","SCAF11")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"DFNA5"]*Coef[1]+dat1[,"GPX4"]*Coef[2]+dat1[,"SCAF11"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6154
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_1_cindex <- data.frame("cohort"='Xiao-Wei Fu',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

#2.Jongmin Kim 
Coef <- c(-0.333,-0.4,0.339,0.387)
gene <- c("CDH1","ID2","MMP9","TCF3")
intersect(gene,colnames(rt)) #都在
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CDH1"]*Coef[1]+dat1[,"ID2"]*Coef[2]+
  dat1[,"MMP9"]*Coef[3]+dat1[,"TCF3"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6489
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_2_cindex <- data.frame("cohort"='Jongmin Kim',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

#3.Zhen Zhang
Coef <- c(-0.13681,-0.07452,-0.09026,-0.1212)
gene <- c("CAT","EHHADH","ALDH5A1","SLC27A5")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CAT"]*Coef[1]+dat1[,"EHHADH"]*Coef[2]+dat1[,"ALDH5A1"]*Coef[3]+
  dat1[,"SLC27A5"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6082
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_3_cindex <- data.frame("cohort"='Zhen Zhang',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

#4.Yue-ling Peng
Coef <- c(0.0029,0.0086,0.0024,0.9751)
gene <- c("HSP90AA1","PPIA","SQSTM1","USP21")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"HSP90AA1"]*Coef[1]+dat1[,"PPIA"]*Coef[2]+dat1[,"SQSTM1"]*Coef[3]+
  dat1[,"USP21"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6061
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_4_cindex <- data.frame("cohort"='Yue-ling Peng',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)

#5.Xiangkun Wang
Coef <- c(-0.5,-0.592,-0.616)
gene <- c("JAK2","STAT5A","STAT6")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"JAK2"]*Coef[1]+dat1[,"STAT5A"]*Coef[2]+
  dat1[,"STAT6"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.657
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_5_cindex <- data.frame("cohort"='Xiangkun Wang',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#6.Shanshan Lu
Coef <- c(0.267389,0.489842)
gene <- c("ANO10", "CLCN2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ANO10"]*Coef[1]+dat1[,"CLCN2"]*Coef[2]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.64
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_6_cindex <- data.frame("cohort"='Shanshan Lu',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#7.Jia-heng Xing
Coef <- c(0.0310975867942454,0.0385177274200257,0.014946834324991)
gene <- c("AKT1", "MAPK3" ,"CASP3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"AKT1"]*Coef[1]+dat1[,"MAPK3"]*Coef[2]+
  dat1[,"CASP3"]*Coef[3]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6557
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_7_cindex <- data.frame("cohort"='Jia-heng Xing',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#8.Hong Peng
Coef<-c(0.0634,0.12544,0.0481,0.0483)
gene <- c("ATP7A", "MTF1", "GLS", "CDKN2A")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ATP7A"]*Coef[1]+dat1[,"MTF1"]*Coef[2]+
  dat1[,"GLS"]*Coef[3]+dat1[,"CDKN2A"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6535
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_8_cindex <- data.frame("cohort"='Hong Peng',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#9.Huaxiang Wang
Coef<-c(0.0618,0.0538,0.1552,0.0428,-0.3199,0.1208)
gene <- c("BUB3", "IGF2BP3","RBM3", "ILF3", "ZC3H13", "CCT3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"BUB3"]*Coef[1]+dat1[,"IGF2BP3"]*Coef[2]+
  dat1[,"RBM3"]*Coef[3]+dat1[,"ILF3"]*Coef[4]+
  dat1[,"ZC3H13"]*Coef[5]+dat1[,"CCT3"]*Coef[6]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6467
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_9_cindex <- data.frame("cohort"='Huaxiang Wang',"C-index"=cindex["C"],
                            "Lower"=lower_cindex,"Upper"=upper_cindex)
#10.Yuqin Tang
Coef<-c(1.3,0.68,1.17,1.49,1.17)
gene <- c("PSRC1", "SOCS2", "TMEM45A", "CCT5", "STC2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"PSRC1"]*Coef[1]+dat1[,"SOCS2"]*Coef[2]+
  dat1[,"TMEM45A"]*Coef[3]+dat1[,"CCT5"]*Coef[4]+
  dat1[,"STC2"]*Coef[5]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6425
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_10_cindex <- data.frame("cohort"='Yuqin Tang',"C-index"=cindex["C"],
                             "Lower"=lower_cindex,"Upper"=upper_cindex)

GSE14520_all_Cindex <- bind_rows(list(GSE14520_MPCDI_cindex,GSE14520_1_cindex,GSE14520_2_cindex
                                  ,GSE14520_3_cindex,GSE14520_4_cindex,GSE14520_5_cindex,
                                  GSE14520_6_cindex,GSE14520_7_cindex,GSE14520_8_cindex,
                                  GSE14520_9_cindex,GSE14520_10_cindex))
GSE14520_all_Cindex <-GSE14520_all_Cindex[order(GSE14520_all_Cindex$C.index,decreasing = T),]
write.csv(GSE14520_all_Cindex,file = "./GSE14520_all_Cindex.csv",row.names = F)
#Plot---------------
data <- read.csv("./GSE14520_all_Cindex.csv",header = T)
data$cohort <- factor(data$cohort,levels = rev(data$cohort))
dat1 <- gather(data,C,value,1:3)
pdf(file = "./GSE14520.pdf",width = 5,height = 4)
ggplot(dat1,aes(cohort,value)) +
  #scale_y_continuous(breaks  = c(0.5,0.6,0.7),limits = c(0.4,0.7))+
  geom_segment(data=data,aes(x=cohort,y=Lower,xend=cohort,yend=Upper),
               color="#175676",size=1,alpha=0.5)+
  geom_point(data=subset(dat1,C=="C.index"),aes(color=C),size=6)+
  scale_color_manual(values = c("steelblue"))+
  geom_hline(aes(yintercept=0.6),linetype=5,col="#175676",alpha=0.5)+
  coord_flip()+
  labs(title = "GSE14520",x="",y="")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none")
dev.off()

#GSE116174
load("~/GSE116174_exp_clin.rda")
rt <- GSE116174_exp_clin[,-c(1:3)]
ALL_Cindex <- read.csv("./All_Cindex.csv",header = T,row.names = 1)
#MPCDI
GSE116174_MPCDI_cindex <- ALL_Cindex[3,]
GSE116174_MPCDI_cindex$cohort <- "MPCDI"
#1.Xiao-Wei Fu  
Coef <- c(0.0182,0.0005,0.0188)
#0.0182*GSDME+0.0005*GPX4+0.0188*SCAF11
gene <- c("DFNA5","GPX4","SFRS2IP")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"DFNA5"]*Coef[1]+dat1[,"GPX4"]*Coef[2]+dat1[,"SFRS2IP"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6154
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_1_cindex <- data.frame("cohort"='Xiao-Wei Fu',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#2.Jongmin Kim 
Coef <- c(-0.333,-0.4,0.339,0.387)
gene <- c("CDH1","ID2","MMP9","TCF3")
intersect(gene,colnames(rt)) #都在
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CDH1"]*Coef[1]+dat1[,"ID2"]*Coef[2]+
  dat1[,"MMP9"]*Coef[3]+dat1[,"TCF3"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6489
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_2_cindex <- data.frame("cohort"='Jongmin Kim',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#3.Zhen Zhang
Coef <- c(-0.13681,-0.07452,-0.09026,-0.1212)
gene <- c("CAT","EHHADH","ALDH5A1","SLC27A5")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CAT"]*Coef[1]+dat1[,"EHHADH"]*Coef[2]+dat1[,"ALDH5A1"]*Coef[3]+
  dat1[,"SLC27A5"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6082
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_3_cindex <- data.frame("cohort"='Zhen Zhang',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#4.Yue-ling Peng
Coef <- c(0.0029,0.0086,0.0024,0.9751)
gene <- c("HSP90AA1","PPIA","SQSTM1","USP21")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"HSP90AA1"]*Coef[1]+dat1[,"PPIA"]*Coef[2]+dat1[,"SQSTM1"]*Coef[3]+
  dat1[,"USP21"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6061
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_4_cindex <- data.frame("cohort"='Yue-ling Peng',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#5.Xiangkun Wang
Coef <- c(-0.5,-0.592,-0.616)
gene <- c("JAK2","STAT5A","STAT6")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"JAK2"]*Coef[1]+dat1[,"STAT5A"]*Coef[2]+
  dat1[,"STAT6"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.657
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_5_cindex <- data.frame("cohort"='Xiangkun Wang',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#6.Shanshan Lu
Coef <- c(0.267389,0.489842)
gene <- c("ANO10", "CLCN2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ANO10"]*Coef[1]+dat1[,"CLCN2"]*Coef[2]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.64
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_6_cindex <- data.frame("cohort"='Shanshan Lu',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#7.Jia-heng Xing
Coef <- c(0.0310975867942454,0.0385177274200257,0.014946834324991)
gene <- c("AKT1", "MAPK3" ,"CASP3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"AKT1"]*Coef[1]+dat1[,"MAPK3"]*Coef[2]+
  dat1[,"CASP3"]*Coef[3]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6557
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_7_cindex <- data.frame("cohort"='Jia-heng Xing',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#8.Hong Peng
Coef<-c(0.0634,0.12544,0.0481,0.0483)
gene <- c("ATP7A", "MTF1", "GLS", "CDKN2A")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ATP7A"]*Coef[1]+dat1[,"MTF1"]*Coef[2]+
  dat1[,"GLS"]*Coef[3]+dat1[,"CDKN2A"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6535
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_8_cindex <- data.frame("cohort"='Hong Peng',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#9.Huaxiang Wang
Coef<-c(0.0618,0.0538,0.1552,0.0428,-0.3199,0.1208)
gene <- c("BUB3", "IGF2BP3","RBM3", "ILF3", "ZC3H13", "CCT3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"BUB3"]*Coef[1]+dat1[,"IGF2BP3"]*Coef[2]+
  dat1[,"RBM3"]*Coef[3]+dat1[,"ILF3"]*Coef[4]+
  dat1[,"ZC3H13"]*Coef[5]+dat1[,"CCT3"]*Coef[6]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6467
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_9_cindex <- data.frame("cohort"='Huaxiang Wang',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#10.Yuqin Tang
Coef<-c(1.3,0.68,1.17,1.49,1.17)
gene <- c("PSRC1", "SOCS2", "TMEM45A", "CCT5", "STC2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"PSRC1"]*Coef[1]+dat1[,"SOCS2"]*Coef[2]+
  dat1[,"TMEM45A"]*Coef[3]+dat1[,"CCT5"]*Coef[4]+
  dat1[,"STC2"]*Coef[5]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6425
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_10_cindex <- data.frame("cohort"='Yuqin Tang',"C-index"=cindex["C"],
                                 "Lower"=lower_cindex,"Upper"=upper_cindex)

GSE116174_all_Cindex <- bind_rows(list(GSE116174_MPCDI_cindex,GSE116174_1_cindex,GSE116174_2_cindex
                                      ,GSE116174_3_cindex,GSE116174_4_cindex,GSE116174_5_cindex,
                                      GSE116174_6_cindex,GSE116174_7_cindex,GSE116174_8_cindex,
                                      GSE116174_9_cindex,GSE116174_10_cindex))
GSE116174_all_Cindex <-GSE116174_all_Cindex[order(GSE116174_all_Cindex$C.index,decreasing = T),]
write.csv(GSE116174_all_Cindex,file = "./GSE116174_all_Cindex.csv",row.names = F)
#Plot---------------
data <- read.csv("./GSE116174_all_Cindex.csv",header = T)
data$cohort <- factor(data$cohort,levels = rev(data$cohort))
dat1 <- gather(data,C,value,1:3)
pdf(file = "./GSE116174.pdf",width = 5,height = 4)
ggplot(dat1,aes(cohort,value)) +
  #scale_y_continuous(breaks  = c(0.5,0.6,0.7),limits = c(0.4,0.7))+
  geom_segment(data=data,aes(x=cohort,y=Lower,xend=cohort,yend=Upper),
               color="#175676",size=1,alpha=0.5)+
  geom_point(data=subset(dat1,C=="C.index"),aes(color=C),size=6)+
  scale_color_manual(values = c("#9370DB"))+
  geom_hline(aes(yintercept=0.6),linetype=5,col="#175676",alpha=0.5)+
  coord_flip()+
  labs(title = "GSE116174",x="",y="")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none")
dev.off()


#ICGC
load("~/ICGC-LIRI-JP/ICGC_exp_clin.rda")
rt <- ICGC_exp_clin[,-c(1,2,5)]
ALL_Cindex <- read.csv("./All_Cindex.csv",header = T,row.names = 1)
#MPCDI
ICGC_MPCDI_cindex <- ALL_Cindex[4,]
ICGC_MPCDI_cindex$cohort <- "MPCDI"
#1.Xiao-Wei Fu  
Coef <- c(0.0182,0.0005,0.0188)
#0.0182*GSDME+0.0005*GPX4+0.0188*SCAF11
gene <- c("DFNA5","GPX4","SCAF11")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"DFNA5"]*Coef[1]+dat1[,"GPX4"]*Coef[2]+dat1[,"SCAF11"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6154
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_1_cindex <- data.frame("cohort"='Xiao-Wei Fu',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#2.Jongmin Kim 
Coef <- c(-0.333,-0.4,0.339,0.387)
gene <- c("CDH1","ID2","MMP9","TCF3")
intersect(gene,colnames(rt)) #都在
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CDH1"]*Coef[1]+dat1[,"ID2"]*Coef[2]+
  dat1[,"MMP9"]*Coef[3]+dat1[,"TCF3"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6489
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_2_cindex <- data.frame("cohort"='Jongmin Kim',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#3.Zhen Zhang
Coef <- c(-0.13681,-0.07452,-0.09026,-0.1212)
gene <- c("CAT","EHHADH","ALDH5A1","SLC27A5")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"CAT"]*Coef[1]+dat1[,"EHHADH"]*Coef[2]+dat1[,"ALDH5A1"]*Coef[3]+
  dat1[,"SLC27A5"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6082
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_3_cindex <- data.frame("cohort"='Zhen Zhang',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#4.Yue-ling Peng
Coef <- c(0.0029,0.0086,0.0024,0.9751)
gene <- c("HSP90AA1","PPIA","SQSTM1","USP21")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"HSP90AA1"]*Coef[1]+dat1[,"PPIA"]*Coef[2]+dat1[,"SQSTM1"]*Coef[3]+
  dat1[,"USP21"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.6061
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_4_cindex <- data.frame("cohort"='Yue-ling Peng',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

#5.Xiangkun Wang
Coef <- c(-0.5,-0.592,-0.616)
gene <- c("JAK2","STAT5A","STAT6")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"JAK2"]*Coef[1]+dat1[,"STAT5A"]*Coef[2]+
  dat1[,"STAT6"]*Coef[3]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance #0.657
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_5_cindex <- data.frame("cohort"='Xiangkun Wang',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#6.Shanshan Lu
Coef <- c(0.267389,0.489842)
gene <- c("ANO10", "CLCN2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ANO10"]*Coef[1]+dat1[,"CLCN2"]*Coef[2]

cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.64
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_6_cindex <- data.frame("cohort"='Shanshan Lu',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#7.Jia-heng Xing
Coef <- c(0.0310975867942454,0.0385177274200257,0.014946834324991)
gene <- c("AKT1", "MAPK3" ,"CASP3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"AKT1"]*Coef[1]+dat1[,"MAPK3"]*Coef[2]+
  dat1[,"CASP3"]*Coef[3]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6557
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_7_cindex <- data.frame("cohort"='Jia-heng Xing',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#8.Hong Peng
Coef<-c(0.0634,0.12544,0.0481,0.0483)
gene <- c("ATP7A", "MTF1", "GLS", "CDKN2A")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"ATP7A"]*Coef[1]+dat1[,"MTF1"]*Coef[2]+
  dat1[,"GLS"]*Coef[3]+dat1[,"CDKN2A"]*Coef[4]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6535
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_8_cindex <- data.frame("cohort"='Hong Peng',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#9.Huaxiang Wang
Coef<-c(0.0618,0.0538,0.1552,0.0428,-0.3199,0.1208)
gene <- c("BUB3", "IGF2BP3","RBM3", "ILF3", "ZC3H13", "CCT3")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"BUB3"]*Coef[1]+dat1[,"IGF2BP3"]*Coef[2]+
  dat1[,"RBM3"]*Coef[3]+dat1[,"ILF3"]*Coef[4]+
  dat1[,"ZC3H13"]*Coef[5]+dat1[,"CCT3"]*Coef[6]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6467
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_9_cindex <- data.frame("cohort"='Huaxiang Wang',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)
#10.Yuqin Tang
Coef<-c(1.3,0.68,1.17,1.49,1.17)
gene <- c("PSRC1", "SOCS2", "TMEM45A", "CCT5", "STC2")
intersect(gene,colnames(rt))
dat1 <- rt[,c("status","time",gene)]
dat1$RS <-dat1[,"PSRC1"]*Coef[1]+dat1[,"SOCS2"]*Coef[2]+
  dat1[,"TMEM45A"]*Coef[3]+dat1[,"CCT5"]*Coef[4]+
  dat1[,"STC2"]*Coef[5]
cindex <- summary(coxph(Surv(time,status)~RS,data = dat1))$concordance#0.6425
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_10_cindex <- data.frame("cohort"='Yuqin Tang',"C-index"=cindex["C"],
                                 "Lower"=lower_cindex,"Upper"=upper_cindex)

ICGC_all_Cindex <- bind_rows(list(ICGC_MPCDI_cindex,ICGC_1_cindex,ICGC_2_cindex
                                      ,ICGC_3_cindex,ICGC_4_cindex,ICGC_5_cindex,
                                      ICGC_6_cindex,ICGC_7_cindex,ICGC_8_cindex,
                                      ICGC_9_cindex,ICGC_10_cindex))
ICGC_all_Cindex <-ICGC_all_Cindex[order(ICGC_all_Cindex$C.index,decreasing = T),]
write.csv(ICGC_all_Cindex,file = "./ICGC_all_Cindex.csv",row.names = F)

#Plot---------------
data <- read.csv("./ICGC_all_Cindex.csv",header = T)
data$cohort <- factor(data$cohort,levels = rev(data$cohort))
dat1 <- gather(data,C,value,1:3)
pdf(file = "./ICGC.pdf",width = 5,height = 4)
ggplot(dat1,aes(cohort,value)) +
  #scale_y_continuous(breaks  = c(0.5,0.6,0.7),limits = c(0.4,0.7))+
  geom_segment(data=data,aes(x=cohort,y=Lower,xend=cohort,yend=Upper),
               color="#175676",size=1,alpha=0.5)+
  geom_point(data=subset(dat1,C=="C.index"),aes(color=C),size=6)+
  scale_color_manual(values = c("#EEAD0E"))+
  geom_hline(aes(yintercept=0.6),linetype=5,col="#175676",alpha=0.5)+
  coord_flip()+
  labs(title = "ICGC",x="",y="")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none")
dev.off()
