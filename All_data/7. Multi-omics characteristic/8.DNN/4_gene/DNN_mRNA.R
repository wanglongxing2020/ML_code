library(tidyverse)
library(pacman)
p_load(limma,edgeR,pheatmap)

load("~/tumordata_374.rda")
data <- tumordata
data[,1:ncol(data)] <- apply(data[,1:ncol(data)],2,as.numeric)
colnames(data) <- substr(colnames(data),1,16)
group <- sapply(strsplit(colnames(data),"\\-"),'[',4)
group <- sapply(strsplit(group,""),'[',1)
group <- gsub("2","1",group)
grouplist <- ifelse(as.numeric(group) == 0,"tumor","normal")
data<-rbind(grouplist,data)
data<-t(data)
data <- as.data.frame(data)
row <- rownames(data)
data[,2:ncol(data)] <- apply(data[,2:ncol(data)],2,as.numeric)

Normaldata <- data[data$`1`=='normal',]
Tumordata <- data[data$`1`=='tumor',]
load("~/TCGA_RS.rda")
all <-TCGA_RS[,-4]
colnames(all)[3] <- "RS"
Tumordata <- Tumordata[substr(rownames(Tumordata),1,12)%in%rownames(all),]
Tumordata <- Tumordata[substr(rownames(Tumordata),15,15)=="1",]

low_sample <- rownames(all)[all$RS=="low"]
high_sample <- rownames(all)[all$RS=="high"]

lowdata <- Tumordata[substr(rownames(Tumordata),1,12)%in%low_sample,]
highdata <- Tumordata[substr(rownames(Tumordata),1,12)%in%high_sample,]

newdata <- rbind(lowdata,highdata)
newdata <- t(newdata)
newdata <- newdata[-1,]
GeneExp <- newdata[,1:ncol(newdata)]
TCGA <- matrix(as.numeric(as.matrix(GeneExp)),nrow = nrow(GeneExp),
               dimnames = list(rownames(GeneExp),colnames(GeneExp)))
TCGA=avereps(TCGA)
names <- colnames(TCGA)
a=as.numeric(substr(names,14,15))
table(a)

connumber =263
treatnumber = 102

expMatrix <- TCGA
colSums(expMatrix)
group_list <- c(rep('Low',263),rep('High',102))
group_list <- factor(group_list,levels = c('Low',"High"),ordered = F)
#表达矩阵校正
exprSet <- expMatrix
boxplot(exprSet,outline=F,notch=T,col=group_list,las=2)
library(limma)
exprSet <-normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=F,notch=T,col=group_list,las=2)
#判断数据是否需要转换
#exprSet <- log2(exprSet+1)
#差异分析
dat <- exprSet
design=model.matrix(~factor(group_list))#High为分子 Low为分母
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef = 2,adjust.method = 'BH')
deg=topTable(fit,coef = 2,adjust.method = 'BH',number = Inf)
deg <- na.omit(deg)

#
deg_gene <- deg[deg$adj.P.Val < 0.05,]
deg_gene <-  deg_gene[order(abs(deg_gene$logFC),decreasing = T),]
write.csv(deg_gene,file = "./DNN_alldeg_mRNA(more).csv") #考虑cmap需不需要

deg_fuji <- deg_gene[abs(deg_gene$logFC) >1,]
write.csv(deg_fuji,file = "./deg_fuji.csv")   #只管前面富集分析，后续DNN分析因为有验证集存在，需要另外取交集

deg_100 <- head(deg_gene,100)
write.csv(deg_100,file = "./DNN_100deg_mRNA(more).csv")


data <- TCGA


deg_100 <- read.csv("./DNN_100deg_mRNA(more).csv",header = T,row.names = 1)

DNN_100mRNA_expr <- data[rownames(deg_100),]
DNN_100mRNA_expr <- t(DNN_100mRNA_expr)
DNN_100mRNA_expr <- as.matrix(DNN_100mRNA_expr)
save(DNN_100mRNA_expr,file = "./DNN_100mRNA_expr.rda")  #得pass看GEO
