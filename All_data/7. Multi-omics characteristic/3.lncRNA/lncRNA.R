library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
load("~/TCGA-LIHC_mRNA.rda")
names(assays(data))
rowdata <- rowData(data)
names(rowdata)
table(rowdata$gene_type)
lnc <-data[rowdata$gene_type == "lncRNA",]
expr_tpm_lnc <- assay(lnc,"tpm_unstrand")
symbol_lnc <- rowData(lnc)$gene_name
head(symbol_lnc)
expr_tpm_lncrna_symbol <- cbind(data.frame(symbol_lnc),
                                 as.data.frame(expr_tpm_lnc))
lnc_expr <- aggregate(.~symbol_lnc,expr_tpm_lncrna_symbol,max)
rownames(lnc_expr) <- lnc_expr$symbol_lnc
lnc_expr <-lnc_expr[,-1]
save(lnc_expr,file = "lnc_expr_tpm.rda")

library(pacman)
p_load(limma,edgeR,pheatmap)
#High-OPCDI
load("~/lnc_expr_tpm.rda")
data <- lnc_expr
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
treatnumber = 102

expMatrix <- TCGA
expMatrix <- log2(expMatrix+1)
colSums(expMatrix)
group_list <- c(rep('Normal',50),rep('Tumor',102))
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
write.csv(deg_gene,file = "./High_alldeg_lncRNA.csv")
deg_gene2 <- deg_gene[abs(deg_gene$logFC) >1,]
write.csv(deg_gene2,file = "./High_FC1deg_lncRNA.csv")


#Low-OPCDI
load("~/lnc_expr_tpm.rda")
data <- lnc_expr
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
treatnumber = 263
#tpm转化为tmp值
expMatrix <- TCGA
expMatrix <- log2(expMatrix+1)
colSums(expMatrix)
group_list <- c(rep('Normal',50),rep('Tumor',263))
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
write.csv(deg_gene,file = "./Low_alldeg_lncRNA.csv")
deg_gene2 <- deg_gene[abs(deg_gene$logFC) >1,]
write.csv(deg_gene2,file = "./Low_FC1deg_lncRNA.csv")
