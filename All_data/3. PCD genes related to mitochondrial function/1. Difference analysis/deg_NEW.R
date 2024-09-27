library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
getGDCprojects()$project_id
LIHC <- GDCquery(project = "TCGA-LIHC", data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "STAR - Counts")
GDCdownload(LIHC,method="api")

data <- GDCprepare(LIHC)
save(data,file = "./TCGA-LIHC_mRNA.rda")
geneexp <- assay(data,i="tpm_unstrand")

rowdata <- rowData(data)
head(rowdata$gene_name)

se_mrna <- data[rowdata$gene_type=="protein_coding",]
se_mrna

expr_tpm_mrna<- assay(se_mrna,"tpm_unstrand")
expr_tpm_mrna[1:10,1:2]

symbol_mrna <- rowData(se_mrna)$gene_name
head(symbol_mrna)

expr_tpm_mrna_symbol <- cbind(data.frame(symbol_mrna),as.data.frame(expr_tpm_mrna))
library(tidyverse)

LIHC_expr <- aggregate(.~symbol_mrna,expr_tpm_mrna_symbol,max)
rownames(LIHC_expr) <- LIHC_expr$symbol_mrna
LIHC_expr <- LIHC_expr[,-1]
save(LIHC_expr,file = "LIHC_expr_tpm.rda")

library(pacman)
p_load(limma,edgeR,pheatmap)

dat <- LIHC_expr
colnames(dat) <- substr(colnames(dat),1,16)

#过滤地表达基因
keep_feature <- rowSums(dat>1) >= 2
table(keep_feature)
dat_filt <- dat[keep_feature,]
dat <- dat_filt

#进行log2(TPM+1)标准化处理
dat <- log2(dat+1)


group <- sapply(strsplit(colnames(dat),"\\-"),'[',4)
group <- sapply(strsplit(group,""),'[',1)
group <- gsub("2","1",group)
grouplist <- ifelse(as.numeric(group) == 0,"tumor","normal")
dat<-rbind(grouplist,dat)
dat<-t(dat)
dat <- as.data.frame(dat)
Normaldat <- dat[dat$`1`=='normal',]
Tumordat <- dat[dat$`1`=='tumor',]

tumordata <- Tumordat[,-1]
tumordata <- t(tumordata) %>%as.data.frame()
save(tumordata,file = "tumordata_374.rda")


#新加
load("~/2.ssGSEA+WGCNA/TCGA_tumor_tpm.rda")
Tumordat <- Tumordat[substr(rownames(Tumordat),1,12)%in%colnames(TCGA),]
Tumordat <- Tumordat[substr(rownames(Tumordat),15,15)==1,]


newdat <- rbind(Normaldat,Tumordat)
newdat <- t(newdat)
newdat <- newdat[-1,]
GeneExp <- newdat[,1:ncol(newdat)]
TCGA <- matrix(as.numeric(as.matrix(GeneExp)),nrow = nrow(GeneExp),
               dimnames = list(rownames(GeneExp),colnames(GeneExp)))
#TCGA=avereps(TCGA) 前面已经去重了
names <- colnames(TCGA)
a=as.numeric(substr(names,14,15))
table(a)

connumber =50
treatnumber = 365

expMatrix <- TCGA
colSums(expMatrix)
group_list <- c(rep('Normal',50),rep('Tumor',365))
group_list <- factor(group_list,levels = c('Normal',"Tumor"),ordered = F)
#limma包默认是排序靠后的 vs 排序靠前的！
#表达矩阵校正
exprSet <-expMatrix 
boxplot(exprSet,outline=F,notch=T,col=group_list,las=2)
library(limma)
exprSet <-normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=F,notch=T,col=group_list,las=2)
#判断数据是否需要转换
#exprSet <- log2(exprSet+1) 不需要
#差异分析
dat <- exprSet
design=model.matrix(~factor(group_list))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef = 2,adjust.method = 'BH')
deg=topTable(fit,coef = 2,adjust.method = 'BH',number = Inf)
deg <- na.omit(deg)

#共16921基因

deg_gene <- deg[deg$adj.P.Val < 0.05,]
write.csv(deg_gene,file = "./LIHC_10995deg.csv")

deg$Change <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 0.5, 
                     ifelse(deg$logFC> 0.5 ,'UP','DOWN'),
                     'NOT')
deg <- deg[deg$Change!="NOT",]

write.csv(deg,file = "./LIHC_4529deg.csv")

table(deg$Change)
#DOWN 2284  UP 2245


library('ggplot2')
p <- ggplot(# 数据、映射、颜色
  deg, aes(x = logFC, y = -log10(adj.P.Val), colour=Change)) +
  geom_point(alpha=0.5, size=3.5) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  #辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  labs(x="log2(Fold change)",y="-log10 (adj.p-value)")+   # 坐标轴# 坐标轴和图标题title="Volcano plot",
  theme_bw()+    #去除背景色
  theme(panel.grid = element_blank())+  #去除网格线
  #xlim(-2, 2)+   #设置坐标轴范围
  #图例
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="top", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
        axis.title.x =element_text(size=18), 
        axis.title.y=element_text(size=18),
        axis.text=element_text(size=14,face = "bold"))

pdf("hotplot_new.pdf", width = 8, height = 8)
p
dev.off()
