library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
#miRNA用的ucsc中下载的log2(RPM+1)格式
#High-OPCDI
files <- list.files(path = "./",
                         pattern = "tsv.gz",
                         full.names = T,
                         recursive = T)
tmp <- readr::read_tsv(files[1])


library(pacman)
p_load(limma,DESeq2,pheatmap)

data <- as.data.frame(tmp)
rownames(data) <- data$miRNA_ID
data <- data[,-1]
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
high <- TCGA_RS[TCGA_RS$RS_group=="high",]
Tumordata <- Tumordata[substr(rownames(Tumordata),1,12)%in%rownames(high),]

newdata <- rbind(Normaldata,Tumordata)
newdata <- t(newdata)
newdata <- newdata[-1,]
GeneExp <- newdata[,1:ncol(newdata)]
TCGA <- matrix(as.numeric(as.matrix(GeneExp)),nrow = nrow(GeneExp),
               dimnames = list(rownames(GeneExp),colnames(GeneExp)))
TCGA=avereps(TCGA)
names <- colnames(TCGA)
a=as.numeric(substr(names,14,15))
table(a)

connumber =50
treatnumber = 101

expMatrix <- TCGA

colSums(expMatrix)
group_list <- c(rep('Normal',50),rep('Tumor',101))
group_list <- factor(group_list,levels = c('Normal',"Tumor"),ordered = F)
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
design=model.matrix(~factor(group_list))#Tumor为分子 Normal为分母
# colnames(design) <- levels(factor(group_list))  #Tumor为分子 Normal为分母
# rownames(design) <- colnames(dat)
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef = 2,adjust.method = 'BH')
deg=topTable(fit,coef = 2,adjust.method = 'BH',number = Inf)
deg <- na.omit(deg)


deg_gene <- deg[deg$adj.P.Val < 0.05,]
deg_gene <-  deg_gene[order(deg_gene$logFC,decreasing = T),]
write.csv(deg_gene,file = "./High_alldeg_miRNA.csv")
deg_gene2 <- deg_gene[abs(deg_gene$logFC) >1,]
write.csv(deg_gene2,file = "./High_FC1deg_miRNA.csv")

#Low-OPCDI
library(pacman)
p_load(limma,DESeq2,pheatmap)

data <- as.data.frame(tmp)
rownames(data) <- data$miRNA_ID
data <- data[,-1]
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
low <- TCGA_RS[TCGA_RS$RS_group=="low",]
Tumordata <- Tumordata[substr(rownames(Tumordata),1,12)%in%rownames(low),]
Tumordata <- Tumordata[substr(rownames(Tumordata),15,15)=="1",]

newdata <- rbind(Normaldata,Tumordata)
newdata <- t(newdata)
newdata <- newdata[-1,]
GeneExp <- newdata[,1:ncol(newdata)]
TCGA <- matrix(as.numeric(as.matrix(GeneExp)),nrow = nrow(GeneExp),
               dimnames = list(rownames(GeneExp),colnames(GeneExp)))
TCGA=avereps(TCGA)
names <- colnames(TCGA)
a=as.numeric(substr(names,14,15))
table(a)

connumber =50
treatnumber = 260
expMatrix <- TCGA

colSums(expMatrix)
group_list <- c(rep('Normal',50),rep('Tumor',260))
group_list <- factor(group_list,levels = c('Normal',"Tumor"),ordered = F)
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
design=model.matrix(~factor(group_list))#Tumor为分子 Normal为分母
# colnames(design) <- levels(factor(group_list))  #Tumor为分子 Normal为分母
# rownames(design) <- colnames(dat)
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef = 2,adjust.method = 'BH')
deg=topTable(fit,coef = 2,adjust.method = 'BH',number = Inf)
deg <- na.omit(deg)

deg_gene <- deg[deg$adj.P.Val < 0.05,]
deg_gene <-  deg_gene[order(deg_gene$logFC,decreasing = T),]
write.csv(deg_gene,file = "./Low_alldeg_miRNA.csv")
deg_gene2 <- deg_gene[abs(deg_gene$logFC) >1,]
write.csv(deg_gene2,file = "./Low_FC1deg_miRNA.csv")
