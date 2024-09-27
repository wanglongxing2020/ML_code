library(GEOquery)
library(limma)
library(umap)
library(tidyverse)
library(stringr)


#GSE14520    488样本--242--221       #用GPL3921
gset <- getGEO("GSE14520", GSEMatrix =TRUE, AnnotGPL=TRUE,destdir = '.',getGPL=T)[[1]]
save(gset,file = "./GSE14520_gset.rda")
gsms <- "undefined"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex+1) }
ex <- as.data.frame(ex)
tmp <- gset@featureData@data
tmp <- tmp[,c("ID","Gene symbol")]
tmp$`Gene symbol` <- gsub("\\///.*","",tmp$`Gene symbol`)
colnames(tmp)[2] <- "Gene_Symbol"
tmp1 <- ex%>%as.data.frame()%>%mutate(ID = rownames(ex))%>%select(ID,everything())
tmp2 <- merge(tmp,tmp1,by="ID")
tmp2 <- tmp2[,-1]

#tmp2 <- tmp2[tmp2$Gene_Symbol%in%deg_gene,]
tmp2 <- tmp2[tmp2$Gene_Symbol!="",]
tmp3 <- aggregate(.~Gene_Symbol,tmp2,max)
GSE14520_exp <- tmp3
GSE14520_exp <- GSE14520_exp[complete.cases(GSE14520_exp),]
save(GSE14520_exp,file = "GSE14520_exp.rda")


#只选肿瘤样本  242例


clin <- read.table("./GSE14520_Extra_Supplement.txt",
                       header = T,na.strings = "NA",sep = "\t",fill =T)
clin <- clin[clin$Tissue.Type == "Tumor",]
clin1 <- clin[,c(3,8,9,15,19,20)]
clin1 <- na.omit(clin1)  #242
clin1$time <- as.numeric(clin1$Survival.months)/12
clin1 <- na.omit(clin1)
clin1$Stage <- gsub("[A-C]","",clin1$TNM.staging)
clin1$status <- clin1$Survival.status

GSE14520_clin <- clin1[,c(1,2,3,8,9,7)]
colnames(GSE14520_clin)[1] <- "Sample"
rownames(GSE14520_clin) <- GSE14520_clin$Sample
GSE14520_clin <- GSE14520_clin[,-1]
save(GSE14520_clin,file = "GSE14520_clin.rda")


rownames(GSE14520_exp) <- GSE14520_exp$Gene_Symbol
GSE14520_exp <- GSE14520_exp[,-1]
exp <- t(GSE14520_exp)
exp_clin <- merge(GSE14520_clin,exp,by="row.names") %>%distinct(Row.names,.keep_all = T)
rownames(exp_clin) <- exp_clin$Row.names
exp_clin <- exp_clin[,-1]
GSE14520_exp_clin <- exp_clin
save(GSE14520_exp_clin,file = "GSE14520_exp_clin.rda")






#GSE116174    64样本
gset <- getGEO("GSE116174", GSEMatrix =TRUE, getGPL=T,destdir = ".")
save(gset,file = "./GSE116174_gset.rda")
if (length(gset) > 1) idx <- grep("GPL13158", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex+1) }
ex <-as.data.frame(ex)
tmp <- gset@featureData@data
tmp <- tmp[,c("ID","Gene Symbol")]
tmp$`Gene Symbol` <- gsub("\\///.*","",tmp$`Gene Symbol`)
colnames(tmp)[2] <- "Gene_Symbol"
tmp1 <- ex%>%as.data.frame()%>%mutate(ID = rownames(ex))%>%select(ID,everything())
tmp2 <- merge(tmp,tmp1,by="ID")
tmp2 <- tmp2[,-1]

#tmp2 <- tmp2[tmp2$Gene_Symbol%in%deg_gene,]
tmp2 <- tmp2[tmp2$Gene_Symbol!="",]
tmp3 <- aggregate(.~Gene_Symbol,tmp2,max)
GSE116174_exp <- tmp3
GSE116174_exp <- GSE116174_exp[complete.cases(GSE116174_exp),]
save(GSE116174_exp,file = "GSE116174_exp.rda")

clin <- read.csv("./GSE116174_HCC64_clin.csv",header = T,row.names = 1,
                 na.strings = "")

clin1 <- clin[,c(2,3,5,9,10)]

ss <- gset@phenoData@data
ss <- ss[,c(2,8)]
rownames(ss) <- ss$source_name_ch1
clin2 <- merge(ss,clin1,by="row.names")
rownames(clin2) <- clin2$geo_accession
clin2 <- clin2[,-c(1:3)]
colnames(clin2) <- c("Gender","Age","Stage","time","status")
clin2$time <- clin2$time/12
clin2$Stage <- ifelse(clin2$Stage == "Ⅰ","I",clin2$Stage)
clin2$Stage[clin2$Stage=="I"] <- "I"
clin2$Stage[clin2$Stage=="Ⅱ"] <- "II"
clin2$Stage[clin2$Stage=="Ⅲ"] <- "III"

GSE116174_clin <- clin2[sort(rownames(clin2)),]
save(GSE116174_clin,file = "GSE116174_clin.rda")

rownames(GSE116174_exp) <- GSE116174_exp$Gene_Symbol
GSE116174_exp <- GSE116174_exp[,-1]
exp <- t(GSE116174_exp)
exp_clin <- merge(GSE116174_clin,exp,by="row.names") %>%distinct(Row.names,.keep_all = T)
rownames(exp_clin) <- exp_clin$Row.names
exp_clin <- exp_clin[,-1]
GSE116174_exp_clin <- exp_clin
save(GSE116174_exp_clin,file = "GSE116174_exp_clin.rda")




#GSE76427    167样本---115
gset <- getGEO("GSE76427", GSEMatrix =TRUE, getGPL=T,destdir = ".")
save(gset,file = "./GSE76427_gset.rda")
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex+1) }
ex <- as.data.frame(ex)
tmp <- gset@featureData@data
tmp <- tmp[,c("ID","Symbol")]
tmp$`Symbol` <- gsub("\\///.*","",tmp$`Symbol`)
colnames(tmp)[2] <- "Gene_Symbol"
tmp1 <- ex%>%as.data.frame()%>%mutate(ID = rownames(ex))%>%select(ID,everything())
tmp2 <- merge(tmp,tmp1,by="ID")
tmp2 <- tmp2[,-1]

#tmp2 <- tmp2[tmp2$Gene_Symbol%in%deg_gene,]
tmp2 <- tmp2[tmp2$Gene_Symbol!="",]
tmp3 <- aggregate(.~Gene_Symbol,tmp2,max)
GSE76427_exp <- tmp3
GSE76427_exp <- GSE76427_exp[complete.cases(GSE76427_exp),]
save(GSE76427_exp,file = "GSE76427_exp.rda")

clin <- gset@phenoData@data
clin <- clin[-str_which(clin$title,"non"),]
ss <- clin[,c(42:52)]
clin1 <-ss[,c(1,4,6,8,11)]
clin1$Age <- clin1$`age (years):ch1`
clin1$Gender <- ifelse(clin1$`gender (1=m, 2=f):ch1`==1,"m","f")
clin1$Stage <- clin1$`tnm_staging_clinical:ch1`
clin1$Stage <- gsub("[A-B]","",clin1$Stage)
clin1$status <- clin1$`event_os:ch1`
clin1$time <- clin1$`duryears_os:ch1`
clin2 <- clin1[,c(6:10)]
GSE76427_clin <- clin2[sort(rownames(clin2)),]
save(GSE76427_clin,file = "GSE76427_clin.rda")


rownames(GSE76427_exp) <- GSE76427_exp$Gene_Symbol
GSE76427_exp <- GSE76427_exp[,-1]
exp <- t(GSE76427_exp)
exp_clin <- merge(GSE76427_clin,exp,by="row.names") %>%distinct(Row.names,.keep_all = T)
rownames(exp_clin) <- exp_clin$Row.names
exp_clin <- exp_clin[,-1]
GSE76427_exp_clin <- exp_clin
save(GSE76427_exp_clin,file = "GSE76427_exp_clin.rda")
