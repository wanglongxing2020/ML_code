library(tidyverse)
library(dplyr)
library(estimate)
library(stringr)

library(Hmisc)
load("~/TCGA_tumor_tpm.rda")
tmp1 <- TCGA

estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(tmp1,file = input.f,sep = '\t',quote = F)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="affymetrix")   ## platform
  scores=read.table(output.ds,skip = 2,header = T,check.names = F)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro="LIHC"
scores  <- estimate(tmp1,pro)
scores <- as.data.frame(scores)
rownames(scores) <- gsub("\\.","-",rownames(scores))
write.csv(scores,file = "LIHC_estimate.csv")

#相关性分析
load("~/TCGA_RS.rda")
scores2 <- scores[rownames(TCGA_RS),]
dat <- TCGA_RS[,4,drop=F]

dat2 <- cbind(dat,scores2)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res<-rcorr(as.matrix(dat2),type="pearson")

#StromalScore:1.195165e-01;ImmuneScore:3.898213e-05;
#ESTIMATEScore:1.081960e-03;TumorPurity:2.084177e-04
res1 <- flattenCorrMatrix(res$r, res$P)
res2 <- res1[res1$row=="RS",]
res2$p_value <- ifelse(res2$p < 0.05,ifelse(res2$p < 0.01,
                                            ifelse(res2$p < 0.001,
                                                   ifelse(res2$p < 0.0001,"****","***"),"**"),"*")
                       ,NA)

write.csv(res2,file = "ESTIMATE_cor_P.csv")

#热图格式
library(pheatmap)
library(RColorBrewer)

dat3 <- dat2[order(dat2$RS,decreasing = F),]
dat3$StromalScore <- 2*((dat3$StromalScore-min(dat3$StromalScore))/(max(dat3$StromalScore)-min(dat3$StromalScore)))-1
dat3$ImmuneScore <- 2*((dat3$ImmuneScore-min(dat3$ImmuneScore))/(max(dat3$ImmuneScore)-min(dat3$ImmuneScore)))-1
dat3$ESTIMATEScore <- 2*((dat3$ESTIMATEScore-min(dat3$ESTIMATEScore))/(max(dat3$ESTIMATEScore)-min(dat3$ESTIMATEScore)))-1
dat3$TumorPurity <- 2*((dat3$TumorPurity-min(dat3$TumorPurity))/(max(dat3$TumorPurity)-min(dat3$TumorPurity)))-1

write.csv(dat3,file = "estimate_hot.csv")


