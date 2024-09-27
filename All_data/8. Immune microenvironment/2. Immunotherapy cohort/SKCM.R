library(GEOquery)
library(limma)
library(umap)
library(tidyverse)
library(survival)
library(survminer)
gset <- getGEO("GSE78220", GSEMatrix =TRUE, getGPL=T,destdir = ".")
clin <-gset[["GSE78220_series_matrix.txt.gz"]]@phenoData
clin <- clin@data


clin1 <- clin[,c(1,2,50:70)]
clin2 <- clin1[,c(1,2,9,16,23)]
colnames(clin2) <- c("title","Sample","Response","time","status")

clin2$Response <- ifelse(clin2$Response=="Progressive Disease","NR","R")
clin2$status <- ifelse(clin2$status=="Alive",0,1)
write.csv(clin2,file = "SKCM_clin.csv",row.names = T)

#基因表达矩阵
clin3 <- clin2[-7,]
data <- read.csv("./GSE78220_FPKM.csv",header = T,row.names = 1)
colnames(data)<- gsub("\\..*$", "",colnames(data))

data2 <- t(data) %>%as.data.frame()
data2 <- data2[-7,]

gene <- read.csv("./gene_4.csv",header = T,row.names = 1)
Gene <- rownames(gene)
data2 <- data2[,Gene]
data2 <- log2(data2+1)

clin4 <- clin3[,c(1,4,5)]
rownames(clin4) <- clin4$title;clin4 <- clin4[,-1]
clin4$time <- clin4$time/365
data2 <- data2[rownames(clin4),]

rt <- cbind(clin4,data2)
rt$RS <-rt$FYN*gene$Coef[1]+rt$HMOX1*gene$Coef[2]+
  rt$LGALS3*gene$Coef[3]+rt$S100A9*gene$Coef[4]
all <- rt

res.cut=surv_cutpoint(all,time = "time",event = "status",variables = c("RS"))
res.cat=surv_categorize(res.cut)
res.cat <-res.cat[clin3$title,]

result <- cbind(clin3,res.cat)

write.csv(result,"./SKCM_result.csv")

table(result$Response,result$RS)
test <- matrix(c(8,3,5,11),nrow = 2,ncol = 2)

chisq.test(test) #p-value = 0.08409

library(ggplot2)
library(ggsci)
res2 <- as.data.frame(table(result$Response,result$RS))
res2$percent <-c(round(8/11,3),round(3/11,3),
                 round(5/16,3),round(11/16,3))
res2$text <- c("72.7%","27.3%","31.2%","68.8%")


pdf("./SKCM.pdf",width = 6,height = 8)
ggplot(res2,aes(x=Var2,y=percent,fill=Var1))+
  geom_bar(stat = 'identity',width = 0.5,colour='black')+
  scale_fill_jama()+
  theme_classic()+xlab("RS")+ylab("Percent weight")+ggtitle("SKCM")+
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  guides(fill=guide_legend(nrow = 1,title = "Response"))+
  geom_text(aes(label=text),position =  position_stack(vjust = 0.5))
dev.off()


#生存分析
library(ggplot2)
library(survival)
library(survminer)

surv <- result
surv <- surv[,-c(4,5)]
fit <- survfit(Surv(time,status) ~ RS,data = surv)
res.cox <- coxph(Surv(time,status) ~RS,data = surv)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS,data = surv)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")
pdf("./SKCM-surv.pdf",width = 6,height = 6)
ggsurvplot(fit,
           data = surv,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =0.5,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama",
           ggtheme = theme_bw(),
           legend="top",
           risk.table = T,
           risk.table.col = 'RS',
           risk.table.y.text = F
) + labs(caption = hr)
dev.off()
