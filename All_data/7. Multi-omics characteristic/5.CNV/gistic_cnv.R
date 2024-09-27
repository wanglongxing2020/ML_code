library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
query <- GDCquery(
  project = "TCGA-LIHC", 
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",
)
GDCdownload(query)
GDCprepare(query, save = T,save.filename = "TCGA-LIHC_CNV.Rdata")
a = load(file = "./TCGA-LIHC_CNV.Rdata")
data1 <- data
tumorCNV  = eval(parse(text = a))
tumorCNV = tumorCNV[,2:7]
tumorCNV = tumorCNV[,c('Sample','Chromosome','Start','End','Num_Probes',"Segment_Mean")]
load("~/TCGA_RS.rda")
rt <- TCGA_RS
rownames(rt) <- paste0(rownames(rt),"-01A")

high_sample = rownames(rt[rt$RS_group=='high',])
low_sample = rownames(rt[rt$RS_group=='low',])
high_CNV = tumorCNV[substring(tumorCNV$Sample,1,16) %in% high_sample,]
low_CNV = tumorCNV[substring(tumorCNV$Sample,1,16) %in% low_sample,]

write.table(high_CNV,file = './High_CNV.txt',sep = '\t',quote = F,row.names = F)
write.table(low_CNV,file = './Low_CNV.txt',sep = '\t',quote = F,row.names = F)

marker_file = read.delim('./snp6.na35.remap.hg38.subset.txt')
marker_file = marker_file[marker_file$freqcnv==FALSE,]
marker_file = marker_file[,1:3]
colnames(marker_file) = c('Marker Name','Chromosome','Marker Position')
write.table(marker_file,file = './marker_file.txt',sep = '\t',quote = F,row.names = F)