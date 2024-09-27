library(tidyverse)
methy_data <- data.table::fread("TCGA-LIHC.methylation450.tsv",data.table = F)
colnames(methy_data)[1] <- "sample"  
rownames(methy_data) <- methy_data[,1]
methy_data <- methy_data[,-1]
table(substr(colnames(methy_data),14,16))


clinical.all <- read.delim("clinical.tsv",header = T,stringsAsFactors = F)
clinical <- clinical.all[,c("case_submitter_id","tissue_or_organ_of_origin")]
table(clinical$tissue_or_organ_of_origin)
origin <- c("Liver")
clinical <- clinical[clinical$tissue_or_organ_of_origin %in% origin,]

sample.all <- read.delim("sample.tsv",header = T,stringsAsFactors = F)
sample_sample <- sample.all[,c("case_submitter_id","sample_submitter_id","sample_type")]
table(sample_sample$sample_type)
tissue_type <- c("Primary Tumor","Solid Tissue Normal")
sample_sample <- sample_sample[sample_sample$sample_type %in% tissue_type,]
merge_info <- merge(sample_sample,clinical,by="case_submitter_id")
merge_info <- merge_info[!duplicated(merge_info$sample_submitter_id),] 
methy_data <- methy_data[,colnames(methy_data) %in% merge_info$sample_submitter_id]
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(methy_data), ]
table(substr(colnames(methy_data),14,16))

#High-OPCDI
load("~/TCGA_RS.rda")
high <- TCGA_RS[TCGA_RS$RS_group=="high",]
high <- paste(rownames(high),"-01A",sep = "")
norm <- colnames(methy_data)[substr(colnames(methy_data),14,16)=="11A"]
Sample <- c(norm,high)
methy_data <- methy_data[,colnames(methy_data) %in% Sample]
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(methy_data), ]

library(ChAMP)
library(dplyr)
library(tibble)
methy_sort <- methy_data[,merge_info$sample_submitter_id]
beta_value = as.matrix(methy_sort)

#beta信号值矩阵里面不能有NA值
library(impute)
beta=impute.knn(beta_value)
sum(is.na(beta))
beta = beta$data
beta = beta+0.00001

#save(beta,file="beta.rda")
#save(merge_info,file="merge.rda")
load("beta.rda")
load("merge.rda")
beta <- as.matrix(beta)
#过滤

myload <- champ.filter(beta = beta,pd = merge_info)
dim(myload$beta)
#save(myload,file = './myload.rda')
#标准化前的QC
QC = champ.QC(beta = myload$beta,pheno = myload$pd$sample_type)
s1 <- table(colnames(myload$beta))
s1 <- as.data.frame(s1)
write.csv(s1,file = "s1.csv")
s1 <- c("TCGA-BC-A10S-11A","TCGA-BC-A110-01A")
keep1 <- setdiff(colnames(myload$beta),s1)
myload$beta <- myload$beta[,keep1]
myload$pd <- myload$pd[-which(myload$pd$sample_submitter_id%in%s1),]
QC = champ.QC(beta = myload$beta,pheno = myload$pd$sample_type)
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(myload$beta),]
#数据标准化（归一化）

myNorm <- champ.norm(beta = myload$beta,arraytype = "450K",cores = 4)
dim(myNorm)
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
table(num.na)
library(stringr)
names(num.na) <- colnames(myNorm)
dt <- names(num.na[num.na > 0])
keep <- setdiff(colnames(myNorm),dt)
myNorm <- myNorm[,keep]
dim(myNorm)
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(myNorm),]

x <- merge_info$sample_type
x[which(x=="Primary Tumor")] <- "Tumor"
x[which(x=="Solid Tissue Normal")] <- "Normal"
merge_info$sample_type <- x
group_list <- merge_info$sample_type

myDMP <- champ.DMP(beta = myNorm,pheno = group_list)
head(myDMP$Tumor_to_Normal)
df_DMP <- myDMP$Tumor_to_Normal

df_DMP <- df_DMP[df_DMP$gene!="",]
DF_0.4_dmp <- df_DMP[which(abs(df_DMP$logFC) >= 0.4),]
DF_0.4_dmp <-DF_0.4_dmp[DF_0.4_dmp$adj.P.Val <0.05, ]
write.csv(DF_0.4_dmp,file ="LIHC_High_DMP.csv")


#最开始再跑一遍
#Low-OPCDI
load("~/TCGA_RS.rda")
low <- TCGA_RS[TCGA_RS$RS_group=="low",]
low <- paste(rownames(low),"-01A",sep = "")
norm <- colnames(methy_data)[substr(colnames(methy_data),14,16)=="11A"]
Sample <- c(norm,low)
methy_data <- methy_data[,colnames(methy_data) %in% Sample]
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(methy_data), ]

library(ChAMP)
library(dplyr)
library(tibble)
methy_sort <- methy_data[,merge_info$sample_submitter_id]
beta_value = as.matrix(methy_sort)

#beta信号值矩阵里面不能有NA值
library(impute)
beta=impute.knn(beta_value)
sum(is.na(beta))
beta = beta$data
beta = beta+0.00001

save(beta,file="beta.rda")
save(merge_info,file="merge.rda")
load("beta.rda")
load("merge.rda")
beta <- as.matrix(beta)
#过滤

myload <- champ.filter(beta = beta,pd = merge_info)
dim(myload$beta)
save(myload,file = './myload.rda')
#标准化前的QC
QC = champ.QC(beta = myload$beta,pheno = myload$pd$sample_type)
s2 <- table(colnames(myload$beta))
s2 <- as.data.frame(s2)
write.csv(s2,file = "s2.csv")
s2 <- c("TCGA-ED-A627-01A","TCGA-CC-5258-01A","TCGA-CC-A7IE-01A",
        "TCGA-UB-A7ME-01A","TCGA-WX-AA47-01A","TCGA-G3-A7M8-01A",
        "TCGA-MR-A520-01A","TCGA-BC-A10X-01A","TCGA-DD-A1ED-01A")
keep1 <- setdiff(colnames(myload$beta),s2)
myload$beta <- myload$beta[,keep1]
myload$pd <- myload$pd[-which(myload$pd$sample_submitter_id%in%s2),]
QC = champ.QC(beta = myload$beta,pheno = myload$pd$sample_type)
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(myload$beta),]

#检查图
#数据标准化（归一化）

myNorm <- champ.norm(beta = myload$beta,arraytype = "450K",cores = 4)
dim(myNorm)
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
table(num.na)
library(stringr)
names(num.na) <- colnames(myNorm)
dt <- names(num.na[num.na > 0])
keep <- setdiff(colnames(myNorm),dt)
myNorm <- myNorm[,keep]
dim(myNorm)
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(myNorm),]

x <- merge_info$sample_type
x[which(x=="Primary Tumor")] <- "Tumor"
x[which(x=="Solid Tissue Normal")] <- "Normal"
merge_info$sample_type <- x
group_list <- merge_info$sample_type

myDMP <- champ.DMP(beta = myNorm,pheno = group_list)
head(myDMP$Tumor_to_Normal)
df_DMP <- myDMP$Tumor_to_Normal
df_DMP <- df_DMP[df_DMP$gene!="",]
DF_0.4_dmp <- df_DMP[which(abs(df_DMP$logFC) >= 0.4),]
DF_0.4_dmp <-DF_0.4_dmp[DF_0.4_dmp$adj.P.Val <0.05, ]
write.csv(DF_0.4_dmp,file ="LIHC_Low_DMP.csv")
