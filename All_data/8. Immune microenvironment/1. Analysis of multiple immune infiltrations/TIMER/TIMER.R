
library(tidyverse)
library(Hmisc)
#数据准备
load("~/TCGA_tumor_tpm.rda")
data <- as.matrix(TCGA)
library(IOBR) 

data[1:4,1:4]

im_timer <- deconvo_tme(eset = data
                        ,method = "timer"
                        ,group_list = rep("lihc",dim(data)[2])
)
write.csv(im_timer,file = "./TIMER_result.csv")

#相关性分析
load("~/TCGA_RS.rda")
res <- as.data.frame(im_timer)
rownames(res) <- res$ID;res<- res[,-1]
res2 <- res[rownames(TCGA_RS),]
dat <- TCGA_RS[,4,drop=F]

dat2 <- cbind(dat,res2)
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

res1 <- flattenCorrMatrix(res$r, res$P)
res2 <- res1[res1$row=="RS",]
res2$p_value <- ifelse(res2$p < 0.05,ifelse(res2$p < 0.01,
                                            ifelse(res2$p < 0.001,
                                                   ifelse(res2$p < 0.0001,"****","***"),"**"),"*")
                       ,NA)
write.csv(res2,file = "TIMER_cor_P.csv")

dat3 <- dat2[order(dat2$RS,decreasing = F),]

write.csv(dat3,file = "TIMER_hot.csv")
