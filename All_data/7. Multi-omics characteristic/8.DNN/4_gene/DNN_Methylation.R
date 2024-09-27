library(tidyverse)
methy_data <- data.table::fread("../methylation/TCGA-LIHC.methylation450.tsv",data.table = F)
colnames(methy_data)[1] <- "sample"  
rownames(methy_data) <- methy_data[,1]
methy_data <- methy_data[,-1]
table(substr(colnames(methy_data),14,16))


clinical.all <- read.delim("../methylation/clinical.tsv",header = T,stringsAsFactors = F)
clinical <- clinical.all[,c("case_submitter_id","tissue_or_organ_of_origin")]
table(clinical$tissue_or_organ_of_origin)
origin <- c("Liver")
clinical <- clinical[clinical$tissue_or_organ_of_origin %in% origin,]

sample.all <- read.delim("../methylation/sample.tsv",header = T,stringsAsFactors = F)
sample_sample <- sample.all[,c("case_submitter_id","sample_submitter_id","sample_type")]
table(sample_sample$sample_type)
tissue_type <- c("Primary Tumor","Solid Tissue Normal")
sample_sample <- sample_sample[sample_sample$sample_type %in% tissue_type,]
merge_info <- merge(sample_sample,clinical,by="case_submitter_id")
merge_info <- merge_info[!duplicated(merge_info$sample_submitter_id),] 
methy_data <- methy_data[,colnames(methy_data) %in% merge_info$sample_submitter_id]
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(methy_data), ]
table(substr(colnames(methy_data),14,16))


load("~/TCGA_RS.rda")
high <- TCGA_RS[TCGA_RS$RS_group=="high",]
high <- paste(rownames(high),"-01A",sep = "")
low <- TCGA_RS[TCGA_RS$RS_group=="low",]
low <- paste(rownames(low),"-01A",sep = "")
Sample <- c(low,high)
methy_data <- methy_data[,colnames(methy_data) %in% Sample]
merge_info <- merge_info[merge_info$sample_submitter_id %in% colnames(methy_data), ]
merge_info$sample_type <- ifelse(merge_info$sample_submitter_id %in% high,"High","Low")



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
s1 <- table(colnames(myload$beta))
s1 <- as.data.frame(s1)
write.csv(s1,file = "s1.csv")
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
x[which(x=="Low")] <- "Low"
x[which(x=="High")] <- "High"
merge_info$sample_type <- x
group_list <- merge_info$sample_type

myDMP <- champ.DMP(beta = myNorm,pheno = group_list)
head(myDMP$Low_to_High)
df_DMP <- myDMP$Low_to_High
df_DMP <- df_DMP[df_DMP$gene!="",]
df_DMP2 <-df_DMP[df_DMP$adj.P.Val < 0.05, ]
write.csv(df_DMP2,file ="DNN_DMP.csv")

deg_gene <-  df_DMP2[order(abs(df_DMP2$logFC),decreasing = T),]
deg_100 <- head(deg_gene,100)
write.csv(deg_100,file = "./DNN_100deg_Methylation.csv")

DNN_100Methylation_expr <- beta[rownames(deg_100),]
DNN_100Methylation_expr <- as.matrix(DNN_100Methylation_expr)
DNN_100Methylation_expr <- DNN_100Methylation_expr-0.00001
table(is.na(DNN_100Methylation_expr))
save(DNN_100Methylation_expr,file = "./DNN_100Methylation_expr.rda")
