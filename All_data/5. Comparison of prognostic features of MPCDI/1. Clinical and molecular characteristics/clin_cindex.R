#TCGA
library(tidyverse)
library(maftools)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(TCGAbiolinks)

clin <- read.csv("./clin_all.csv",header = T,row.names = 1)
colnames(clin)[1] <- "Tumor_Sample_Barcode"

query_SNV <- GDCquery(project = "TCGA-LIHC",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query_SNV)
mafFilePath = dir(path = "./GDCdata/TCGA-LIHC/",pattern = "masked.maf.gz$",full.names = T,recursive=T)
x = lapply(mafFilePath, data.table::fread, skip = "Hugo_Symbol")
x = data.table::rbindlist(l = x, use.names = TRUE, fill = TRUE)
lam1 = maftools::read.maf(maf = x,clinicalData = clin,isTCGA = T)
##计算TMB
tmb <- tmb(maf = lam1,
           captureSize = 35,  #tcga使用的是GRch38参考基因组，长度约35Mb
           logScale = T)
tmb <- tmb[,c(1,3)]
colnames(tmb)[2] <- "TMB"
LIHC_tmb <- tmb
write.csv(LIHC_tmb,file = "LIHC_tmb.csv")

#TP53、KRAS
#提取数据
snpdata <-lam1@data
tidydata <- snpdata[,c("Hugo_Symbol","Variant_Type",
                       "Variant_Classification","Tumor_Sample_Barcode")]
#删除无义突变
tidydata <- tidydata[which(tidydata$Variant_Classification != "Nonsense_Mutation"),
                     c("Hugo_Symbol","Variant_Type","Tumor_Sample_Barcode")]
library(reshape2)
rsdata = dcast(tidydata,Hugo_Symbol~Tumor_Sample_Barcode)
rownames(rsdata)<-rsdata$Hugo_Symbol
rsdata <- rsdata[,-1]#0代表野生型，其他数字代表有突变
rsdata1 <- rsdata[c("TP53","KRAS"),]
rsdata2 <- t(rsdata1)%>%as.data.frame()
rsdata2$Sample <- rownames(rsdata2)
LIHC_2mutation <- rsdata2
write.csv(LIHC_2mutation,file = "LIHC_2mutation.csv")

load("~/TCGA_RS.rda")
clin <- read.csv("./clin_all.csv",header = T,row.names = 1)
TCGA_RS$Sample <- rownames(TCGA_RS)
TCGA_RS <- TCGA_RS[,c("Sample","RS","RS_group")]
rownames(TCGA_RS) <- NULL
TCGA_clin <- merge(clin,TCGA_RS,by="Sample")
TCGA_clin$time <- TCGA_clin$time/365
save(TCGA_clin,file = "./TCGA_clin.rda")


ALL_Cindex <- read.csv("./All_Cindex.csv",header = T,row.names = 1)
#MPCDI
TCGA_MPCDI_cindex <- ALL_Cindex[1,]
TCGA_MPCDI_cindex$cohort <- "MPCDI"

#TMB 
tmb <- read.csv("./LIHC_tmb.csv",header = T,row.names = 1)
colnames(tmb)[1] <- "Sample" 
data1 <- merge(TCGA_clin,tmb,by="Sample")
cindex <- summary(coxph(Surv(time,status)~TMB,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_TMB_cindex <- data.frame("cohort"='TMB',"C-index"=cindex["C"],
                              "Lower"=lower_cindex,"Upper"=upper_cindex)


#TP53、KRAS
rsdata2 <- read.csv("LIHC_2mutation.csv",header = T,row.names = 1)
data2 <- merge(rsdata2,TCGA_clin,by="Sample")

cindex <- summary(coxph(Surv(time,status)~TP53,data = data2))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_TP53_cindex <- data.frame("cohort"='TP53',"C-index"=cindex["C"],
                               "Lower"=lower_cindex,"Upper"=upper_cindex)
cindex <- summary(coxph(Surv(time,status)~KRAS,data = data2))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_KRAS_cindex <- data.frame("cohort"='KRAS',"C-index"=cindex["C"],
                               "Lower"=lower_cindex,"Upper"=upper_cindex)

#MSI
#来自lihc_tcga_pan_can_atlas_2018_clinical_data
MSI <-read.csv("./LIHC_MSI.csv",header = T)


data3 <- merge(MSI,TCGA_clin,by="Sample")
cindex <- summary(coxph(Surv(time,status)~MSI,data = data3))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_MSI_cindex <- data.frame("cohort"='MSI',"C-index"=cindex["C"],
                              "Lower"=lower_cindex,"Upper"=upper_cindex)

#Age 
data1 <- TCGA_clin 
cindex <- summary(coxph(Surv(time,status)~Age,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_Age_cindex <- data.frame("cohort"='Age',"C-index"=cindex["C"],
                              "Lower"=lower_cindex,"Upper"=upper_cindex)

#Gender 
data1 <- TCGA_clin
cindex <- summary(coxph(Surv(time,status)~Gender,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_Gender_cindex <- data.frame("cohort"='Gender',"C-index"=cindex["C"],
                                 "Lower"=lower_cindex,"Upper"=upper_cindex)

#T 
data1 <- TCGA_clin
data_T <- data1[data1$T!='TX',]
data_T <- data_T[data_T$T!="'--",]
cindex <- summary(coxph(Surv(time,status)~T,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_T_cindex <- data.frame("cohort"='T',"C-index"=cindex["C"],
                             "Lower"=lower_cindex,"Upper"=upper_cindex)

#M 
data1 <- TCGA_clin
data_M <- data1[data1$M!='MX',]
cindex <- summary(coxph(Surv(time,status)~M,data = data_M))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_M_cindex <- data.frame("cohort"='M',"C-index"=cindex["C"],
                             "Lower"=lower_cindex,"Upper"=upper_cindex)

#N 
data1 <- TCGA_clin
data_N <- data1[data1$N!='NX',]
data_N <- data_N[data_N$N != "'--",]
cindex <- summary(coxph(Surv(time,status)~N,data = data_N))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_N_cindex <- data.frame("cohort"='N',"C-index"=cindex["C"],
                             "Lower"=lower_cindex,"Upper"=upper_cindex)


#Stage 
data1 <- TCGA_clin
data1 <- data1[data1$Stage!="'--",]
cindex <- summary(coxph(Surv(time,status)~Stage,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_Stage_cindex <- data.frame("cohort"='Stage',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

TCGA_clin_Cindex <- bind_rows(list(TCGA_MPCDI_cindex,TCGA_Age_cindex,TCGA_Gender_cindex
                                   ,TCGA_T_cindex,TCGA_M_cindex,TCGA_N_cindex,
                                   TCGA_Stage_cindex,TCGA_TMB_cindex,TCGA_MSI_cindex,
                                   TCGA_TP53_cindex,TCGA_KRAS_cindex))
write.csv(TCGA_clin_Cindex,file = "./TCGA_clin_Cindex.csv",row.names = F)

#GSE14520
load("~/GSE14520_RS.rda")
load("~/GSE14520_clin.rda")

#MPCDI
GSE14520_MPCDI_cindex <- ALL_Cindex[2,]
GSE14520_MPCDI_cindex$cohort <- "MPCDI"

GSE14520_RS <- GSE14520_RS[,c("RS","RS_group")]
GSE14520_clin <- merge(GSE14520_clin,GSE14520_RS,by="row.names")
save(GSE14520_clin,file = "./GSE14520_clin.rda")

#Age 
data1 <- GSE14520_clin
cindex <- summary(coxph(Surv(time,status)~Age,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_Age_cindex <- data.frame("cohort"='Age',"C-index"=cindex["C"],
                              "Lower"=lower_cindex,"Upper"=upper_cindex)


#Gender 
data1 <- GSE14520_clin
cindex <- summary(coxph(Surv(time,status)~Gender,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_Gender_cindex <- data.frame("cohort"='Gender',"C-index"=cindex["C"],
                                 "Lower"=lower_cindex,"Upper"=upper_cindex)

#Stage 
data1 <- GSE14520_clin
data1 <- data1[data1$Stage!=".",]
cindex <- summary(coxph(Surv(time,status)~Stage,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE14520_Stage_cindex <- data.frame("cohort"='Stage',"C-index"=cindex["C"],
                                "Lower"=lower_cindex,"Upper"=upper_cindex)

GSE14520_clin_Cindex <- bind_rows(list(GSE14520_MPCDI_cindex,GSE14520_Age_cindex,
                                       GSE14520_Gender_cindex,
                                       GSE14520_Stage_cindex))
write.csv(GSE14520_clin_Cindex,file = "./GSE14520_clin_Cindex.csv",row.names = F)


#GSE116174
load("~/GSE116174_RS.rda")
load("~/GSE116174_clin.rda")

#MPCDI
GSE116174_MPCDI_cindex <- ALL_Cindex[3,]
GSE116174_MPCDI_cindex$cohort <- "MPCDI"

GSE116174_RS <- GSE116174_RS[,c("RS","RS_group")]
GSE116174_clin <- merge(GSE116174_clin,GSE116174_RS,by="row.names")
save(GSE116174_clin,file = "./GSE116174_clin.rda")

#Age 
data1 <- GSE116174_clin
cindex <- summary(coxph(Surv(time,status)~Age,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_Age_cindex <- data.frame("cohort"='Age',"C-index"=cindex["C"],
                                  "Lower"=lower_cindex,"Upper"=upper_cindex)

#Gender 
data1 <- GSE116174_clin
cindex <- summary(coxph(Surv(time,status)~Gender,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_Gender_cindex <- data.frame("cohort"='Gender',"C-index"=cindex["C"],
                                     "Lower"=lower_cindex,"Upper"=upper_cindex)

#Stage 
data1 <- GSE116174_clin
cindex <- summary(coxph(Surv(time,status)~Stage,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GSE116174_Stage_cindex <- data.frame("cohort"='Stage',"C-index"=cindex["C"],
                                    "Lower"=lower_cindex,"Upper"=upper_cindex)

GSE116174_clin_Cindex <- bind_rows(list(GSE116174_MPCDI_cindex,GSE116174_Age_cindex,
                                       GSE116174_Gender_cindex,
                                       GSE116174_Stage_cindex))
write.csv(GSE116174_clin_Cindex,file = "./GSE116174_clin_Cindex.csv",row.names = F)


#ICGC
load("~/ICGC_RS.rda")
load("~/ICGC_clin.rda")

#MPCDI
ICGC_MPCDI_cindex <- ALL_Cindex[4,]
ICGC_MPCDI_cindex$cohort <- "MPCDI"

ICGC_RS <- ICGC_RS[,c("RS","RS_group")]
rownames(ICGC_clin) <- ICGC_clin$icgc_donor_id
ICGC_clin <- ICGC_clin[,-1]
ICGC_clin <- merge(ICGC_clin,ICGC_RS,by="row.names")
save(ICGC_clin,file = "./ICGC_clin.rda")

#Age 
data1 <- ICGC_clin
cindex <- summary(coxph(Surv(time,status)~Age,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_Age_cindex <- data.frame("cohort"='Age',"C-index"=cindex["C"],
                                   "Lower"=lower_cindex,"Upper"=upper_cindex)

#Gender 
data1 <- ICGC_clin
cindex <- summary(coxph(Surv(time,status)~Gender,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_Gender_cindex <- data.frame("cohort"='Gender',"C-index"=cindex["C"],
                                      "Lower"=lower_cindex,"Upper"=upper_cindex)


#Stage 
data1 <- ICGC_clin
cindex <- summary(coxph(Surv(time,status)~Stage,data = data1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
ICGC_Stage_cindex <- data.frame("cohort"='Stage',"C-index"=cindex["C"],
                                     "Lower"=lower_cindex,"Upper"=upper_cindex)


ICGC_clin_Cindex <- bind_rows(list(ICGC_MPCDI_cindex,ICGC_Age_cindex,
                                        ICGC_Gender_cindex,
                                        ICGC_Stage_cindex))
write.csv(ICGC_clin_Cindex,file = "./ICGC_clin_Cindex.csv",row.names = F)

#-----------plot------

#TCGA-LUAD
plt1 <- read.csv("./TCGA_clin_Cindex.csv",header = T)
plt1$cohort <- factor(plt1$cohort,levels = plt1$cohort)

pdf("./TCGA.pdf",width = 4,height = 4)
fill_cols <- c(	
  MPCDI="#DC143C",Age="#87CEFA",Gender="#20B2AA",
               T="#4682B4",M="#F4A460",N="#778899",Stage="#66CDAA",TMB="#BDB76B",
               MSI="#FFD700",TP53="#8A2BE2",KRAS="#FF1493")
ggplot(plt1,aes(cohort,C.index))+
  geom_col(width = .5,aes(fill=cohort))+
  scale_fill_manual(values = fill_cols)+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),position=position_dodge(0.9),width=0.15)+
  labs(title = "TCGA-LIHC",x="",y="C-index")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  guides(fill = F)+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.8))+
  scale_x_discrete(expand = c(0.05,0))
dev.off()


#GSE14520
plt2 <- read.csv("./GSE14520_clin_Cindex.csv",header = T)
plt2$cohort <- factor(plt2$cohort,levels = c("MPCDI","Age","Gender","Stage"))
pdf("./GSE14520.pdf",width = 2,height = 4)
fill_cols <- c(MPCDI="#DC143C",Age="#87CEFA",
               Gender="#4682B4",Stage="#778899")
ggplot(plt2,aes(cohort,C.index))+
  geom_col(width = .4,aes(fill=cohort))+
  scale_fill_manual(values = fill_cols)+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.15)+
  labs(title = "GSE14520",x="",y="C-index")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  guides(fill = F)+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.8))
dev.off()

#GSE116174
plt3 <- read.csv("./GSE116174_clin_Cindex.csv",header = T)
plt3$cohort <- factor(plt3$cohort,levels = c("MPCDI","Age","Gender","Stage"))
pdf("./GSE116174.pdf",width = 2,height = 4)
fill_cols <- c(MPCDI="#DC143C",Age="#87CEFA",
               Gender="#4682B4",Stage="#778899")
ggplot(plt3,aes(cohort,C.index))+
  geom_col(width = .4,aes(fill=cohort))+
  scale_fill_manual(values = fill_cols)+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.15)+
  labs(title = "GSE116174",x="",y="C-index")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  guides(fill = F)+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.8))
dev.off()

#ICGC
plt4 <- read.csv("./ICGC_clin_Cindex.csv",header = T)
plt4$cohort <- factor(plt4$cohort,levels = c("MPCDI","Age","Gender","Stage"))
pdf("./ICGC.pdf",width = 2,height = 4)
fill_cols <- c(MPCDI="#DC143C",Age="#87CEFA",
               Gender="#4682B4",Stage="#778899")
ggplot(plt4,aes(cohort,C.index))+
  geom_col(width = .4,aes(fill=cohort))+
  scale_fill_manual(values = fill_cols)+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.15)+
  labs(title = "ICGC",x="",y="C-index")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                   plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  guides(fill = F)+
  scale_y_continuous(expand = c(0,0.01),limits = c(0,0.8))
dev.off()

