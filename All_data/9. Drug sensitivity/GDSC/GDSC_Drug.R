rm(list = ls())
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(stringr)
library(psych)
dir='../DataFiles/Training Data/'
dir(dir)
exp= readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
exp[1:4,1:4]
dim(exp) #17419个基因，805个细胞系


drug = readRDS(file = file.path(dir,"GDSC2_Res.rds")) #805细胞系，198个药物
dim(drug)
#排除含NA值超20%的化合物
#计算每列NA的数量
na_counts <- colSums(is.na(drug))
# 计算每列的总数
total_counts <- nrow(drug)

# 计算每列NA的比例
na_percentage <- na_counts / total_counts

# 找出NA比例超过20%的列
columns_to_remove <- na_percentage > 0.2

# 删除这些列
drug2 <- drug[, !columns_to_remove]
dim(drug2) #805细胞系，179个药物
#采用最近邻(k-NN)方法对缺失值进行估计
library(devtools)
install_github("ltorgo/DMwR2",ref="master")
library(DMwR2)
drug3 <- as.data.frame(drug2)
drug3 <- knnImputation(drug3,k=10)

#重要说明：这里使用了 e^IC50，因为 IC50 已经是实际的 ln 值/log 转换，而保罗#has的calcPhenotype 函数将进行幂转换（我认为最好不要同时进行两种转换）。
drug4 <- exp(as.matrix(drug3))
drug4[1:4,1:4]
identical(rownames(drug2),colnames(exp))

test <- read.table("~/TCGA_LIHC_tpm.txt",header=T, sep="\t", row.names=1,check.names=F)

test <- as.matrix(test)


calcPhenotype(trainingExprData = exp,
              trainingPtype = drug4,
              testExprData = test,
              batchCorrect = 'standardize',  #   "eb" for array,standardize  for rnaseq
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

testPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', row.names = 1,check.names = F)

#差异分析
library(tidyverse)
library(rstatix)
drug <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', row.names = 1,check.names = F)
drug <- t(drug) %>%as.data.frame()

drug <- rownames_to_column(drug,var = "Drug")

load("~/TCGA_RS.rda")
clin <- TCGA_RS[,3,drop=F]
clin$Sample <- rownames(clin)
colnames(clin)[1] <- "group"
clin$group <- ifelse(clin$group=="high","High","Low")
clin1 <- clin[,c("Sample","group")]
#差异分析

df <- drug %>%
  as_tibble()%>%
  pivot_longer(-1,names_to = "Sample",values_to = "AUC")%>%
  left_join(clin1,by=c("Sample" = "Sample"))


#计算FC
dfFC = df %>%
  group_by(Drug,group) %>%  
  summarise(mean = mean(AUC,na.rm=T)) %>%               
  pivot_wider(names_from = group,values_from = mean) %>%  
  summarise(FC = High/Low)                         
dfFC
#计算p值
dfP = df %>%
  group_by(Drug) %>%
  wilcox_test(AUC ~ group) 
dfP

#对P值FDR校正
dfP_FDR = dfP %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  
dfP_FDR

#合并
dfData <- drug %>%
  as_tibble() %>%
  left_join(dfFC) %>%
  left_join(dfP_FDR)
save(dfData,file = "./drugdata.rda")

dfData2 <- dfData[,-c(2:366)]
dfData2$logFC <- log2(dfData2$FC)

#这里我选择的是log2FC<-0.3(high组的AUC值应该选更低的), 筛选得到8个药物
dfData2 <- dfData2[order(dfData2$logFC,decreasing = F),]
dfData3 <-dfData2[dfData2$FDR<0.05,]
df_drug <-dfData3$Drug[dfData3$logFC < -0.3]



#相关分析
clin2 <- TCGA_RS[,c("RS"),drop=F]
drug2 <- drug
rownames(drug2) <- drug2$Drug
drug2 <- drug2[,-1]
drug2 <- t(drug2) %>%as.data.frame()

df2 <- merge(clin2,drug2,by="row.names")
rownames(df2) <- df2$Row.names;df2<-df2[,-1]
df3 <- apply(df2,2,as.numeric)

n<-data.frame()
for (i in 2:ncol(df3)){
  Drug<-colnames(df3)[i]
  cor_value <- cor.test(df3[,1],df3[,i],method = "spearman")
  co <- cor_value$estimate
  names(co) <- NULL
  p <- cor_value$p.value 
  s <- data.frame(Drug=Drug,r=co,p.value=p)
  n <- rbind(s,n)
}
write.csv(n,file = "./cor_result.csv")
cor_result <- n[n$p.value<0.05,]
cor_result <- cor_result[cor_result$r < -0.3,]
cor_drug <- cor_result$Drug  #8个

#最终得到5种药物
inter_drug <- intersect(df_drug,cor_drug)
write.csv(inter_drug,file = "./GDSC_inter_drug.csv")

#作图
library(ggplot2)
library(ggpubr)
library(ggsci)
dat1 <- df[df$Drug%in%inter_drug,]
dat1$Drug <-gsub("_.*$", "", dat1$Drug)

a <- dfData3[dfData3$Drug%in%inter_drug,]
b <- a[match(sort(a$Drug),a$Drug),]


pdf(file = "./GDSC_boxplot.pdf",height = 6,width = 12)
ggboxplot(dat1, x = "Drug", y = "AUC",
          fill = "group")+
  scale_fill_jama()+
  stat_compare_means(aes(group = group),
                     method = 'wilcox.test',
                     size = 3,label = "p.signif",hide.ns = T)+
  xlab("")+ylab("Estimated score")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),axis.ticks.y = element_blank())
dev.off()

dat2 <- cor_result[cor_result$Drug%in%inter_drug,]
dat2$Drug <-gsub("_.*$", "", dat2$Drug)
dat2 <- dat2[order(dat2$r),]
dat2$Drug <- factor(dat2$Drug,levels = rev(dat2$Drug))
dat2$P  <- -log10(dat2$p.value) 


pdf("./GDSC_cor.pdf",width = 10,height = 6)
ggplot(dat2,aes(r,Drug))+
  geom_segment(aes(x = 0,xend=r,y=Drug,yend=Drug),size=0.5,linetype = "dashed")+
  geom_point(shape=21,aes(size=P),fill='#4682B4',color='#4682B4')+
  scale_x_reverse()+
  theme_classic()+
  xlab("Correlation coefficient")+ylab("")+
  theme(legend.position = "bottom",
        axis.line.y = element_blank())+
  guides(size=guide_legend(nrow = 1))+
  labs(size="-log10(P-value)")
dev.off()
