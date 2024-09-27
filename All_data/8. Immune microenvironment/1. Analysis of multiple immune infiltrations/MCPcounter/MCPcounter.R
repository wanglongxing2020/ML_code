install.packages(c("devtools","curl"))
library(devtools)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(tidyverse)
library(Hmisc)
#输入数据
load("~/TCGA_tumor_tpm.rda")
index=order(rowMeans(TCGA),decreasing = T)
TPM_ordered=TCGA[index,]
TPM <- as.matrix(TPM_ordered)

library(MCPcounter)
probesets <- read.table('./probesets.txt',
                        sep="\t",
                        stringsAsFactors=FALSE,
                        colClasses="character")
head(probesets,n=3)
genes <- read.table('./genes.txt',
                    sep="\t",
                    stringsAsFactors=FALSE,
                    header=TRUE,
                    colClasses="character",
                    check.names=FALSE)
head(genes,n=3)

res<- MCPcounter.estimate(TPM,
                          featuresType="HUGO_symbols",
                          probesets=probesets,
                          genes=genes)
head(res[,1:4],n=3)

res <- t(res) %>%as.data.frame()
#相关性分析
load("~/TCGA_RS.rda")
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
write.csv(res2,file = "MCP_cor_P.csv")

dat3 <- dat2[order(dat2$RS,decreasing = F),]

write.csv(dat3,file = "MCP_hot.csv")
