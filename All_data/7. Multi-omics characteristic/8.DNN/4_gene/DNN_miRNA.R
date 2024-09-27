library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)


files <- list.files(path = "~/7. Multi-omics characteristic/miRNA",
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

connumber =260
treatnumber = 101

expMatrix <- TCGA
colSums(expMatrix)
group_list <- c(rep('Low',260),rep('High',101))
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
write.csv(deg_gene,file = "./DNN_alldeg_miRNA.csv")


deg_100 <- head(deg_gene,100)
write.csv(deg_100,file = "./DNN_100deg_miRNA.csv")



data <- TCGA
DNN_100miRNA_expr <- data[rownames(deg_100),]
DNN_100miRNA_expr <- t(DNN_100miRNA_expr)
DNN_100miRNA_expr <- as.matrix(DNN_100miRNA_expr)
save(DNN_100miRNA_expr,file = "./DNN_100miRNA_expr.rda")


