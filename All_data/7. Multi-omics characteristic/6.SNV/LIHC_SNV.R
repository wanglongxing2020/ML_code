library(tidyverse)

load("~/TCGA_RS.rda")
clin<-TCGA_RS
clin$Tumor_Sample_Barcode <- rownames(clin)
clin <- clin[,c(5,3)]
clin <- clin[order(clin$RS_group,decreasing = F),]
colnames(clin)[2] <- "RS"
clin$RS <- ifelse(clin$RS=="high","High","Low")

library(dplyr)
library(maftools)
library(tidyverse)
library(TCGAbiolinks)
query_SNV <- GDCquery(project = "TCGA-LIHC",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query_SNV)
mafFilePath = dir(path = "./GDCdata/TCGA-LIHC/",pattern = "masked.maf.gz$",full.names = T,recursive=T)
x = lapply(mafFilePath, data.table::fread, skip = "Hugo_Symbol")
x = data.table::rbindlist(l = x, use.names = TRUE, fill = TRUE)
x2 <- x
x2$Tumor_Sample_Barcode <- substr(x2$Tumor_Sample_Barcode,1,12)
x2 <-x2[x2$Tumor_Sample_Barcode%in%clin$Tumor_Sample_Barcode,]
laml = maftools::read.maf(maf = x2,clinicalData = clin,isTCGA = T)

#提取颜色
library(ggsci)
pal_jama(palette = c("default"), alpha = 1)(7)
colors <- list(RS=c("High"="#374E55","Low"="#DF8F44"))

save(laml,file = "./laml.rda")

pdf("./oncoplot.pdf",width = 10,height = 8)
oncoplot(maf = laml, top = 20,clinicalFeatures = "RS",sortByAnnotation = T,annotationColor = colors)
dev.off()

pdf("./summary.pdf",width = 10,height = 8)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

High<- subsetMaf(maf = laml,clinQuery="RS%in%'High'",mafObj = T,isTCGA = T)
Low<- subsetMaf(maf = laml,clinQuery="RS%in%'Low'",mafObj = T,isTCGA = T)
# 使用mafCompare函数比较差异突变基因(Fisher精确检验)
comparison <- mafCompare(m1 = Low, m2 = High,m1Name = "Low",m2Name = "High")
res <- comparison[["results"]]
write.csv(res,file = "./Deg_SNVgene.csv")


#TP53: High < Low


# 绘制比较结果
# comparison[["results"]] <- comparison[["results"]][1:10,]
forestPlot(mafCompareRes=comparison,pVal=0.05,color=c("maroon","royalblue"),geneFontSize=0.8)
