library(tidyverse)
# 加载需要使用的R包
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)
#TCGA
load("~/TCGA_exp_clin39.rda")
gene <- read.csv("./TCGA_outTab_16.csv",header = T)
gene <- gene$id
rt <-rt[,c("time","status",gene)]
colnames(rt)[1:2]<- c("OS.time","OS")
rt_new <- rt%>%rownames_to_column(var = "ID")
Training_TCGA <- rt_new
#GSE14520 221样本 都在   
load("~/GSE14520_exp_clin.rda")
rt2 <- GSE14520_exp_clin
rt2 <- rt2[,c("time","status",gene)]
colnames(rt2)[1:2]<- c("OS.time","OS")
rt2_new <- rt2%>%rownames_to_column(var = "ID")
Testing_GSE14520 <- rt2_new
#GSE116174  64样本  CXCL8不在但别名IL8在   "CTSL"不在但别名"CTSL1"在
gene2 <- gene
gene2[4] <- "CTSL1"
gene2[5] <- "IL8"
load("~/GSE116174_exp_clin.rda")
rt3 <- GSE116174_exp_clin
rt3 <- rt3[,c("time","status",gene2)]
colnames(rt3)[6] <- "CTSL"
colnames(rt3)[7] <- "CXCL8"
colnames(rt3)[1:2]<- c("OS.time","OS")
rt3_new <- rt3%>%rownames_to_column(var = "ID")
Testing_GSE116174 <- rt3_new
#ICGC     231样本  "CTSL"  "CXCL8"不在    CTSL1在 IL8在
load("~/ICGC-LIRI-JP/ICGC_exp_clin.rda")
rt4 <- ICGC_exp_clin
rt4 <- rt4[,c("time","status",gene2)]
colnames(rt4)[6] <- "CTSL"
colnames(rt4)[7] <- "CXCL8"
colnames(rt4)[1:2]<- c("OS.time","OS")
rt4_new <- rt4%>%rownames_to_column(var = "ID")
Testing_ICGC <- rt4_new

mm <- list(Training_TCGA = Training_TCGA, Testing_GSE14520 = Testing_GSE14520,
           Testing_GSE116174 = Testing_GSE116174,Testing_ICGC = Testing_ICGC)
save(mm,file = "model_datasets.rda")

result <- data.frame()
est_data <- mm$Training_TCGA
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})
rm(mm)


seed <- 120
set.seed(seed)
rf_nodesize <-5

#1.RSF ####
#1.1 RSF 
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result,cc)

rid <-var.select(fit, verbose = F)$topvars
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
#1.2 RSF + CoxBoost
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF',' + CoxBoost')
result <- rbind(result,cc)
#1.3 RSF + Enet [alpha=0-1] alpha = 1为Lasso;alpha = 0为Ridge
x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF',' + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}

#1.4 RSF + GBM
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF',' + GBM')
result <- rbind(result,cc)
#1.5 RSF + plsRcox
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF',' + plsRcox')
result <- rbind(result,cc)
#1.6 RSF + StepCox [backward/forward/both]
set.seed(seed)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF',' + StepCox','[',direction,']')
  result <- rbind(result,cc)
}
#1.7 RSF + SuperPC
set.seed(seed)
data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF',' + SuperPC')
result <- rbind(result,cc)
#1.8 RSF + survival-SVM
set.seed(seed)
fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)

rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF',' + survival-SVM')
result <- rbind(result,cc)
#####
#2.Enet ####
#2.1 Enet [alpha=0-1] alpha = 1为Lasso;alpha = 0为Ridge
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}
#####
#3.StepCox ####
#3.1.1 StepCox [forward]
for (direction in c("forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}
#3.2.1 StepCox [both]
for (direction in c("both")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}

rid <- names(coef(fit))
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

#3.2.2 StepCox [both] + CoxBoost
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [both]',' + CoxBoost')
result <- rbind(result,cc)
#3.2.3 StepCox [both] + Enet [alpha=0-1] alpha = 1为Lasso;alpha = 0为Ridge
x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox [both]',' + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}
#3.2.4 StepCox [both] + GBM
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [both]',' + GBM')
result <- rbind(result,cc)
#3.2.5 StepCox [both] + plsRcox
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [both]',' + plsRcox')
result <- rbind(result,cc)
#3.2.6 StepCox [both] + RSF
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [both]',' + RSF')
result <- rbind(result,cc)
#3.2.7 StepCox [both] + SuperPC
set.seed(seed)
data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [both]',' + SuperPC')
result <- rbind(result,cc)
#3.2.8 StepCox [both] + survival-SVM
set.seed(seed)
fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)

rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [both]',' + survival-SVM')
result <- rbind(result,cc)

#3.3.1 StepCox [backward]
for (direction in c("backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}

rid <- names(coef(fit))
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

#3.3.2 StepCox [backward] + CoxBoost
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [backward]',' + CoxBoost')
result <- rbind(result,cc)
#3.3.3 StepCox [backward] + Enet [alpha=0-1] alpha = 1为Lasso;alpha = 0为Ridge
x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox [backward]',' + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}
#3.3.4 StepCox [backward] + GBM
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [backward]',' + GBM')
result <- rbind(result,cc)
#3.3.5 StepCox [backward] + plsRcox
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [backward]',' + plsRcox')
result <- rbind(result,cc)
#3.3.6 StepCox [backward] + RSF
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [backward]',' + RSF')
result <- rbind(result,cc)
#3.3.7 StepCox [backward] + SuperPC
set.seed(seed)
data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=2,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [backward]',' + SuperPC')
result <- rbind(result,cc)
#3.3.8 StepCox [backward] + survival-SVM
set.seed(seed)
fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)

rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('StepCox [backward]',' + survival-SVM')
result <- rbind(result,cc)
#####
#4.CoxBoost ####
#4.1 CoxBoost
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)

rid <- names(coef(fit)[abs(coef(fit))>0])
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

#4.2 CoxBoost + RSF
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost',' + RSF')
result <- rbind(result,cc)
#4.3 CoxBoost + Enet [alpha=0-1] alpha = 1为Lasso;alpha = 0为Ridge
x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost',' + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}

#4.4 CoxBoost + GBM
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost',' + GBM')
result <- rbind(result,cc)
#4.5 CoxBoost + plsRcox
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost',' + plsRcox')
result <- rbind(result,cc)
#4.6 CoxBoost + StepCox [backward/forward/both]
set.seed(seed)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost',' + StepCox','[',direction,']')
  result <- rbind(result,cc)
}
#4.7 CoxBoost + SuperPC
set.seed(seed)
data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost',' + SuperPC')
result <- rbind(result,cc)
#4.8 CoxBoost + survival-SVM
set.seed(seed)
fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)

rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost',' + survival-SVM')
result <- rbind(result,cc)
#####
#5.plsRcox ####
#5.1 plsRcox
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,pre_var],time=est_dd$OS.time,status=est_dd$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)
#####
#6.SuperPC ####
#6.1 SuperPC
data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,featurenames=colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)
#####
#7.GBM ####
#7.1 GBM
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)
#####
#8.survival-SVM ####
#8.1 survival-SVM
fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd, gamma.mu = 1)

rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')
result <- rbind(result,cc)
#####
#9.Lasso Enet [alpha=1] ####
#9.1 Lasso Enet已计算
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

set.seed(seed)
alpha=1
fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
fit2 = glmnet(x1,x2,family = "cox", alpha = alpha, lambda = fit$lambda.min)
co = coef(fit2)[, 1]
rid <- names(co)[co!=0]
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

#9.2 Lasso + CoxBoost
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso',' + CoxBoost')
result <- rbind(result,cc)
#9.3 Lasso + RSF
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso',' + RSF')
result <- rbind(result,cc)

#9.4 Lasso + GBM
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso',' + GBM')
result <- rbind(result,cc)
#9.5 Lasso + plsRcox
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso',' + plsRcox')
result <- rbind(result,cc)
#9.6 Lasso + StepCox [backward/forward/both]
set.seed(seed)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso',' + StepCox','[',direction,']')
  result <- rbind(result,cc)
}
#9.7 Lasso + SuperPC
set.seed(seed)
data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso',' + SuperPC')
result <- rbind(result,cc)
#9.8 Lasso + survival-SVM
set.seed(seed)
fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)

rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso',' + survival-SVM')
result <- rbind(result,cc)
#####
result2 <- result
result2$Model <- gsub('α','alpha',result2$Model)
result2$Model <- gsub("-| ", "", result2$Model)
result2$Model <- gsub('Enet\\[alpha=0\\]','Ridge',result2$Model)
result2$Model <- gsub('Enet\\[alpha=1\\]','Lasso',result2$Model)
final <- result2 %>% pivot_wider(id_cols =Model,names_from = ID,values_from = Cindex )
final <- as.data.frame(final)
colnames(final) <- c("Model","TCGA-LIHC","GSE14520","GSE116174","ICGC")
write.csv(final,file = "./Cindex_16_all.csv",row.names = F)

#画图部分 ####
Cindex_mat <- final
rownames(Cindex_mat) <- Cindex_mat$Model
Cindex_mat <- Cindex_mat[,-1]
avg_Cindex <- apply(Cindex_mat, 1, mean) # 计算每种算法在所有队列中平均C-index
avg_Cindex <- sort(avg_Cindex, decreasing = T)     # 对各算法C-index由高到低排序
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      # 对C-index矩阵排序
avg_Cindex <- avg_Cindex[rownames(Cindex_mat)]
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)


CohortCol <- c("firebrick","steelblue","#9370DB","#EEAD0E") # 你可以替换这两种颜色为你喜欢的颜色

names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)


cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 黄绿配色
              #col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
              rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
              cluster_columns = FALSE, cluster_rows = FALSE, # 不进行聚类，无意义
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf("./Cindex_16_new.pdf", width = (cellwidth * ncol(Cindex_mat) + 3) * 2, height = cellheight * nrow(Cindex_mat) * 0.45)

draw(hm)
invisible(dev.off())
#####

#最佳模型 StepCox [both] + RSF

##3.2.6 StepCox [both] + RSF
set.seed(seed)
for (direction in c("both")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}

rid <- names(coef(fit))
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
var.select(fit, verbose = F)$topvars #4个 "S100A9" "FYN"    "LGALS3" "HMOX1"
pdf("./RSF.pdf",width = 10,height = 5)
plot(fit)
dev.off()

#系数获取
Coeffit <- coxph(Surv(time=OS.time,event=OS) ~ .,data = est_dd2)
coef(Coeffit)
coef <- data.frame(Gene=names(coef(Coeffit)),Coef=unname(coef(Coeffit)))
write.csv(coef,file = "./gene_4.csv",row.names = F)
#保存在gene_4.csv
library(ggplot2)
library(tidyverse)
data <- read.csv("./gene_4.csv",header = T)
data$direaction <- ifelse(data$Coef > 0,"P","N")
data$direaction <- as.factor(data$direaction)
pdf("./Coefficient.pdf",width = 10,height = 5)
ggplot(data,aes(Coef,Gene))+geom_segment(aes(x = 0,xend=Coef,y=Gene,yend=Gene),color="black",size=1,alpha=0.5)+
  geom_point(shape=21,size=5,aes(fill=direaction))+geom_vline(xintercept = 0,linetype="dotted",color="black",linewidth=1,alpha=0.5)+
  scale_fill_manual(values = c("N"="#4682B4","P"="#8B0000"))+theme_bw()+ggtitle("")+
  xlab("Coefficient")+ylab("")+theme(plot.background = element_blank(),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank())
dev.off()
#TCGA   
new_data1 <- rt[,c('OS.time','OS',rid)]
new_data1$RS <-new_data1$FYN*coef$Coef[1]+new_data1$HMOX1*coef$Coef[2]+
  new_data1$LGALS3*coef$Coef[3]+new_data1$S100A9*coef$Coef[4]
tcga2 <- new_data1
colnames(tcga2)[1:2] <- c("time","status")
  
cindex <- summary(coxph(Surv(time,status)~RS,data = tcga2))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
TCGA_cindex <- data.frame("cohort"='TCGA',"C-index"=cindex["C"],
                          "Lower"=lower_cindex,"Upper"=upper_cindex)
#GSE14520
new_data2 <- rt2[,c('OS.time','OS',rid)]
new_data2$RS <-new_data2$FYN*coef$Coef[1]+new_data2$HMOX1*coef$Coef[2]+
  new_data2$LGALS3*coef$Coef[3]+new_data2$S100A9*coef$Coef[4]
GEO1 <- new_data2
colnames(GEO1)[1:2] <- c("time","status")

cindex <- summary(coxph(Surv(time,status)~RS,data = GEO1))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GEO1_cindex <- data.frame("cohort"='GSE14520',"C-index"=cindex["C"],
                          "Lower"=lower_cindex,"Upper"=upper_cindex)

#GSE116174
new_data3 <- rt3[,c('OS.time','OS',rid)]
new_data3$RS <-new_data3$FYN*coef$Coef[1]+new_data3$HMOX1*coef$Coef[2]+
  new_data3$LGALS3*coef$Coef[3]+new_data3$S100A9*coef$Coef[4]
GEO2 <- new_data3
colnames(GEO2)[1:2] <- c("time","status")

cindex <- summary(coxph(Surv(time,status)~RS,data = GEO2))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
GEO2_cindex <- data.frame("cohort"='GSE116174',"C-index"=cindex["C"],
                          "Lower"=lower_cindex,"Upper"=upper_cindex)

#ICGC
new_data4 <- rt4[,c('OS.time','OS',rid)]
new_data4$RS <-new_data4$FYN*coef$Coef[1]+new_data4$HMOX1*coef$Coef[2]+
  new_data4$LGALS3*coef$Coef[3]+new_data4$S100A9*coef$Coef[4]
icgc <- new_data4
colnames(icgc)[1:2] <- c("time","status")

cindex <- summary(coxph(Surv(time,status)~RS,data = icgc))$concordance
upper_cindex <- cindex["C"] + 1.96*cindex["se(C)"]
lower_cindex <- cindex["C"] - 1.96*cindex["se(C)"]
icgc_cindex <- data.frame("cohort"='ICGC',"C-index"=cindex["C"],
                          "Lower"=lower_cindex,"Upper"=upper_cindex)

All_Cindex <- bind_rows(list(TCGA_cindex,GEO1_cindex,GEO2_cindex,icgc_cindex))

write.csv(All_Cindex,file = "All_Cindex.csv",row.names = F)

#1,3,5years AUC_ROC plot
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)
library(ggsci)
#TCGA
TCGA <- tcga2[,c("time","status","RS")]


ROC_rt=timeROC(T=TCGA$time,delta = TCGA$status,
               marker = TCGA$RS,cause = 1,
               weighting = "aalen",
               times = c(1,3,5),ROC = TRUE)
pdf(file="TCGA-ROC.pdf",width = 5.5,height = 5.5)
plot(ROC_rt,time=1,col="#FA8072",title=F,lwd=2)
plot(ROC_rt,time=3,col="#FFD700",add=T,title=F,lwd=2)
plot(ROC_rt,time=5,col="#32CD32",add=T,title=F,lwd=2)
title(main="TCGA-LIHC")
legend("bottomright",
       c(paste0("AUC at 1 years: ",sprintf("%.3f",ROC_rt$AUC[1])),
         paste0("AUC at 3 years: ",sprintf("%.3f",ROC_rt$AUC[2])),
         paste0("AUC at 5 years: ",sprintf("%.3f",ROC_rt$AUC[3]))),
       col=c("#FA8072","#FFD700","#32CD32"),lwd = 2,bty = "n")
dev.off()

#GSE14520

GSE14520 <- GEO1[,c("time","status","RS")]
ROC_rt=timeROC(T=GSE14520$time,delta = GSE14520$status,
               marker = GSE14520$RS,cause = 1,
               weighting = "aalen",
               times = c(1,3,5),ROC = TRUE)
pdf(file="GSE14520-ROC.pdf",width = 5.5,height = 5.5)
plot(ROC_rt,time=1,col="#FA8072",title=F,lwd=2)
plot(ROC_rt,time=3,col="#FFD700",add=T,title=F,lwd=2)
plot(ROC_rt,time=5,col="#32CD32",add=T,title=F,lwd=2)
title(main="GSE14520")
legend("bottomright",
       c(paste0("AUC at 1 years: ",sprintf("%.3f",ROC_rt$AUC[1])),
         paste0("AUC at 3 years: ",sprintf("%.3f",ROC_rt$AUC[2])),
         paste0("AUC at 5 years: ",sprintf("%.3f",ROC_rt$AUC[3]))),
       col=c("#FA8072","#FFD700","#32CD32"),lwd = 2,bty = "n")
dev.off()

#GSE116174

GSE116174 <- GEO2[,c("time","status","RS")]
ROC_rt=timeROC(T=GSE116174$time,delta = GSE116174$status,
               marker = GSE116174$RS,cause = 1,
               weighting = "aalen",
               times = c(1,3,5),ROC = TRUE)
pdf(file="GSE116174-ROC.pdf",width = 5.5,height = 5.5)
plot(ROC_rt,time=1,col="#FA8072",title=F,lwd=2)
plot(ROC_rt,time=3,col="#FFD700",add=T,title=F,lwd=2)
plot(ROC_rt,time=5,col="#32CD32",add=T,title=F,lwd=2)
title(main="GSE116174")
legend("bottomright",
       c(paste0("AUC at 1 years: ",sprintf("%.3f",ROC_rt$AUC[1])),
         paste0("AUC at 3 years: ",sprintf("%.3f",ROC_rt$AUC[2])),
         paste0("AUC at 5 years: ",sprintf("%.3f",ROC_rt$AUC[3]))),
       col=c("#FA8072","#FFD700","#32CD32"),lwd = 2,bty = "n")
dev.off()

#ICGC

ICGC <- icgc[,c("time","status","RS")]

ROC_rt=timeROC(T=ICGC$time,delta = ICGC$status,
               marker = ICGC$RS,cause = 1,
               weighting = "aalen",
               times = c(1,3),ROC = TRUE)
pdf(file="ICGC-ROC.pdf",width = 5.5,height = 5.5)
plot(ROC_rt,time=1,col="#FA8072",title=F,lwd=2)
plot(ROC_rt,time=3,col="#FFD700",add=T,title=F,lwd=2)
title(main="ICGC")
legend("bottomright",
       c(paste0("AUC at 1 years: ",sprintf("%.3f",ROC_rt$AUC[1])),
         paste0("AUC at 3 years: ",sprintf("%.3f",ROC_rt$AUC[2]))),
       col=c("#FA8072","#FFD700"),lwd = 2,bty = "n")
dev.off()

#高低CDI组生存分析
library(patchwork)
#TCGA

TCGA <- tcga2[,c(1,2,ncol(tcga2))]

surv <- TCGA
surv$status <-ifelse(surv$status=="0",0,1)
res.cut=surv_cutpoint(surv,time = "time",event = "status",variables = c("RS"))
res.cat=surv_categorize(res.cut)
surv<- res.cat
#这里可保存RS的分组信息
Dat1 <- TCGA[,3,drop=F]
TCGA_RS <- cbind(surv,Dat1)
colnames(TCGA_RS)[3] <- "RS_group"
save(TCGA_RS,file = "TCGA_RS.rda")

fit <- survfit(Surv(time,status) ~ RS,data = surv)
res.cox <- coxph(Surv(time,status) ~RS,data = surv)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS,data = surv)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")
pdf("./TCGA-surv.pdf",width = 5.5,height = 5.5)
ggsurvplot(fit,
           data = surv,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =1,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama") + labs(caption = hr)+
  ggtitle("TCGA-LIHC")

dev.off()

#GSE14520

GSE14520 <- GEO1[,c(1,2,ncol(GEO1))]

surv <- GSE14520
surv$status <-ifelse(surv$status=="0",0,1)
res.cut=surv_cutpoint(surv,time = "time",event = "status",variables = c("RS"))
res.cat=surv_categorize(res.cut)
surv<- res.cat
#
Dat2 <- GSE14520[,3,drop=F]
GSE14520_RS <- cbind(surv,Dat2)
colnames(GSE14520_RS)[3] <- "RS_group"
save(GSE14520_RS,file = "GSE14520_RS.rda")

fit <- survfit(Surv(time,status) ~ RS,data = surv)
res.cox <- coxph(Surv(time,status) ~RS,data = surv)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS,data = surv)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")
pdf("./GSE14520-surv.pdf",width = 5.5,height = 5.5)
ggsurvplot(fit,
           data = surv,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =1,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama") + labs(caption = hr)+
  ggtitle("GSE14520")
dev.off()

#GSE116174
GSE116174 <- GEO2[,c(1,2,ncol(GEO2))]

surv <- GSE116174
surv$status <-ifelse(surv$status=="0",0,1)
res.cut=surv_cutpoint(surv,time = "time",event = "status",variables = c("RS"))
res.cat=surv_categorize(res.cut)
surv<- res.cat
#
Dat3 <- GSE116174[,3,drop=F]
GSE116174_RS <- cbind(surv,Dat3)
colnames(GSE116174_RS)[3] <- "RS_group"
save(GSE116174_RS,file = "GSE116174_RS.rda")


fit <- survfit(Surv(time,status) ~ RS,data = surv)
res.cox <- coxph(Surv(time,status) ~RS,data = surv)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS,data = surv)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")
pdf("./GSE116174-surv.pdf",width = 5.5,height = 5.5)
ggsurvplot(fit,
           data = surv,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =1,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama") + labs(caption = hr)+
  ggtitle("GSE116174")
dev.off()

#ICGC
ICGC <- icgc[,c(1,2,ncol(icgc))]

surv <- ICGC
surv$status <-ifelse(surv$status=="0",0,1)
res.cut=surv_cutpoint(surv,time = "time",event = "status",variables = c("RS"))
res.cat=surv_categorize(res.cut)
surv<- res.cat
#
Dat4 <- ICGC[,3,drop=F]
ICGC_RS <- cbind(surv,Dat4)
colnames(ICGC_RS)[3] <- "RS_group"
save(ICGC_RS,file = "ICGC_RS.rda")

fit <- survfit(Surv(time,status) ~ RS,data = surv)
res.cox <- coxph(Surv(time,status) ~RS,data = surv)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS,data = surv)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")
pdf("./ICGC-surv.pdf",width = 5.5,height = 5.5)
ggsurvplot(fit,
           data = surv,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =1,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama") + labs(caption = hr)+
  ggtitle("ICGC")
dev.off()
