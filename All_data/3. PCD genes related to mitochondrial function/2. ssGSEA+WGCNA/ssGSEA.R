library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(GSEABase)
library(rtracklayer)
gene <- read.csv("DAG.marker.csv",header = T)
gene <- gene$X
gset <- c("Hallmark_my_geneset","NA",gene)
gset <- gset%>% 
  as.data.frame() %>% 
  t()
write.table(gset,file = "Hallmark_my_geneset.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

mygeneset <- "Hallmark_my_geneset.gmt"
# 获取gmt文件内容
geneset <- getGmt(mygeneset)

#处理下表达矩阵
load("~/1.Difference analysis/tumordata_374.rda")
TCGA <- t(tumordata)
TCGA <- as.data.frame(TCGA)
TCGA <- TCGA[substr(rownames(TCGA),15,15)== "1",] #选择“1”的原发肿瘤
TCGA <- t(TCGA) %>%as.data.frame()


#发现有重复的
clin_all <- read.csv("./clin_all.csv",row.names = 1)
colnames(TCGA) <- substr(colnames(TCGA),1,12)
TCGA <- TCGA[,colnames(TCGA)%in% clin_all$Sample]
rowname <- rownames(TCGA)
TCGA <- apply(TCGA,2,as.numeric)
TCGA <- data.frame(TCGA,row.names = rowname)
colnames(TCGA) <- gsub("\\.","-",colnames(TCGA))
save(TCGA,file="TCGA_tumor_tpm.rda")


gsva_matrix<- gsva(as.matrix(TCGA),geneset,method='ssgsea',
                   kcdf='Gaussian',ssgsea.norm=TRUE,abs.ranking=TRUE)


mydata <- as.data.frame(gsva_matrix)
mydata <- as.data.frame(t(mydata))
rownames(mydata) <- gsub("\\.","-",rownames(mydata))
save(mydata,file = "ssGSEA.rda")
write.csv(mydata,file = "ssGSEA.csv")

