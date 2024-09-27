library(tidyverse)
library(Hmisc)


load("~/TCGA_tumor_tpm.rda")
tmp1 <- TCGA
write.table(tmp1,file = "TCGA_LIHC_tpm.txt",sep = "\t",quote = F)
source("Cibersort.R")
#TCGA-LIHC
sig_matrix <- read.table(file = "./CIBERSORT_ref.txt",header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
sig_matrix <- log2(sig_matrix+1) %>%as.data.frame()
write.table(sig_matrix,file = "ref.txt",sep = "\t",quote = F)

TCGA_TME.results <- CIBERSORT(sig_matrix = "ref.txt",mixture_file ="TCGA_LIHC_tpm.txt" ,perm = 1000,QN= T)
write.csv(TCGA_TME.results,file = "CIBERSORT_TME_results.csv")

ciber <- TCGA_TME.results[,1:22]
#相关性分析
load("~/TCGA_RS.rda")
ciber2 <- ciber[rownames(TCGA_RS),]
dat <- TCGA_RS[,4,drop=F]

dat2 <- cbind(dat,ciber2)
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

write.csv(res2,file = "CIBERSORT_cor_P.csv")

dat3 <- dat2[order(dat2$RS,decreasing = F),]

write.csv(dat3,file = "cibersort_hot.csv")
