library(tidyverse)
library(survival)
library(survminer)
library(forestplot)
library(MASS)
library(regplot)
library(rms)
library(timeROC)

load("~/TCGA_clin.rda")

data <- TCGA_clin[,c(1,2,3,4,5,8,11)]

data <- data[data$Stage!="'--",]  #365-----341
colnames(data)[7] <- "RS"
data$Stage <- ifelse(data$Stage =="Stage I",1,ifelse(data$Stage == "Stage II",2,ifelse(data$Stage=="Stage III",3,4)))
data$RS <- ifelse(data$RS == "low",1,2)
uniresult=data.frame()
for (i in colnames(data)[4:7]) {
  unicox <- coxph(Surv(time,status) ~ data[,i],data = data)
  unisum = summary(unicox)
  coxP=unisum$coefficients[,"Pr(>|z|)"]
  
  uniresult=rbind(uniresult,
                  cbind(id=i,
                        HR=round(unisum$conf.int[,"exp(coef)"],2),
                        HR.95L=round(unisum$conf.int[,"lower .95"],2),
                        HR.95H=round(unisum$conf.int[,"upper .95"],2),
                        pvalue=coxP)
  )
}
write.csv(uniresult,"uniresult.csv",row.names = F)
result <- uniresult
result$Type <- c("NS","NS","Risky","Risky")
result[,2:5] <- apply(result[,2:5],2,as.numeric)
result$id[4] <- "MPCDI"
result$id <- factor(result$id,levels = c("MPCDI","Stage","Gender","Age"))
result$'p.val' <- round(result$pvalue,4)
result$'p.val' <- c("0.1253","0.2209","<0.001","<0.001")
result$'logP' <- -log10(result$pvalue)
result$'HR(95%CI)' <-paste0(uniresult$HR,"[",uniresult$HR.95L," to ",uniresult$HR.95H,"]")

col <- c("#808080","#990000")
annotation1 <- data.frame(
  x=rep(c(-0.5,-2,-3.5),each=4),
  y=rep(c(1,2,3,4),times=3),
  label=c('<0.001', '<0.001', '0.2209','0.1253',
          '2.98[2.06 to 4.33]','1.64[1.34 to 2.01]','0.79[0.55 to 1.15]','1.01[1 to 1.03]',
          'MPCDI','Stage','Gender','Age'))

pdf("./unicox.pdf",width = 11,height = 4)
ggplot(result,aes(HR,id))+geom_point(aes(color=Type,size=logP))+
  geom_errorbarh(aes(xmax=HR.95H,xmin=HR.95L),height=0.2)+
  scale_x_continuous(limits = c(-4,5),breaks = c(0,1,3,5),expand = c(0,0))+
  scale_color_manual(values = col)+labs(size = "-log10(p.val)")+
  geom_vline(xintercept = c(0,1))+xlab("HR(95%CI)")+ylab("")+theme_bw()+
  theme(axis.text.y  = element_blank(),axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_rect(aes(xmin=-4,xmax=0,ymin=0,ymax=0.5),fill="#F5DEB3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=0.5,ymax=1.5),fill="#D3D3D3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=1.5,ymax=2.5),fill="#F5DEB3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=2.5,ymax=3.5),fill="#D3D3D3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=3.5,ymax=4.5),fill="#F5DEB3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=4.5,ymax=5.5),fill="#D3D3D3",alpha=0.1)+
  geom_text(data=annotation1,aes(x=x,y=y,label=label))+
  geom_text(x=-0.5,y=5,label="p.val")+
  geom_text(x=-2,y=5,label="HR(95%CI)")
dev.off()

#多因素cox
data1 <- data
data1$Gender <- ifelse(data1$Gender=="female",1,2)
COX <- coxph(Surv(time,status)~Age+Gender+Stage+RS,data = data1)

##筛选Stage+RS
multiresult=data.frame()
multisum = summary(COX)
coxP=multisum$coefficients[,"Pr(>|z|)"]
multiresult=rbind(multiresult,
                  cbind(
                    HR=round(multisum$conf.int[,"exp(coef)"],2),
                    HR.95L=round(multisum$conf.int[,"lower .95"],2),
                    HR.95H=round(multisum$conf.int[,"upper .95"],2),
                    pvalue=coxP)
)

write.csv(multiresult,"multiresult.csv")
result <- multiresult
result$id <- rownames(result)
result <- result%>%dplyr::select(id,everything())
result$Type <- c("NS","NS","Risky","Risky")
result[,2:5] <- apply(result[,2:5],2,as.numeric)
result$id[4] <- "MPCDI"
result$id <- factor(result$id,levels = c("MPCDI","Stage","Gender","Age"))
result$'p.val' <- round(result$pvalue,4)
result$'p.val' <- c("0.3322","0.6634","<0.001","<0.001")
result$'logP' <- -log10(result$pvalue)
result$'HR(95%CI)' <-paste0(multiresult$HR,"[",multiresult$HR.95L," to ",multiresult$HR.95H,"]")

col <- c("#808080","#990000")
annotation1 <- data.frame(
  x=rep(c(-0.5,-2,-3.5),each=4),
  y=rep(c(1,2,3,4),times=3),
  label=c('<0.001', '<0.001', '0.6634','0.3322',
          '2.61[1.79 to 3.83]','1.52[1.24 to 1.86]','0.92[0.63 to 1.34]','1.01[0.99 to 1.02]',
          'MPCDI','Stage','Gender','Age'))
pdf("./multicox.pdf",width = 11,height = 4)
ggplot(result,aes(HR,id))+geom_point(aes(color=Type,size=logP))+
  geom_errorbarh(aes(xmax=HR.95H,xmin=HR.95L),height=0.2)+
  scale_x_continuous(limits = c(-4,5),breaks = c(0,1,3,5),expand = c(0,0))+
  scale_color_manual(values = col)+labs(size = "-log10(p.val)")+
  geom_vline(xintercept = c(0,1))+xlab("HR(95%CI)")+ylab("")+theme_bw()+
  theme(axis.text.y  = element_blank(),axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_rect(aes(xmin=-4,xmax=0,ymin=0,ymax=0.5),fill="#F5DEB3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=0.5,ymax=1.5),fill="#D3D3D3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=1.5,ymax=2.5),fill="#F5DEB3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=2.5,ymax=3.5),fill="#D3D3D3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=3.5,ymax=4.5),fill="#F5DEB3",alpha=0.1)+
  geom_rect(aes(xmin=-4,xmax=0,ymin=4.5,ymax=5.5),fill="#D3D3D3",alpha=0.1)+
  geom_text(data=annotation1,aes(x=x,y=y,label=label))+
  geom_text(x=-0.5,y=5,label="p.val")+
  geom_text(x=-2,y=5,label="HR(95%CI)")
dev.off()



#nomogram
#绘制列线图
dat1 <- data
dat1$Stage <- ifelse(dat1$Stage == 1,"I",
                     ifelse(dat1$Stage == 2,"II",
                            ifelse(dat1$Stage==3,"III","IV")))
dat1$Stage <- factor(dat1$Stage,levels = c("I","II","III","IV"))
dat1$RS <- ifelse(dat1$RS == 1,"Low","High")
dat1$RS <- factor(dat1$RS,levels = c("High","Low"))
colnames(dat1)[7] <- "MPCDI"
res.cox=coxph(Surv(time,status) ~Stage+MPCDI,data=dat1)
sum.surv<-summary(res.cox)

pdf(file="nomogram.pdf",width = 12,height = 5)
regplot(res.cox,
        plots = c("bean","boxes"),
        col.grid=gray(c(0.8,0.95)),
        clickable=F,
        title="Nomogram",
        points=T,
        droplines=T,
        #observation=dat1[1,],  #设置观察者
        rank="sd",
        failtime = c(1,3,5),
        prfail = F,
        dencol="#F5DEB3",
        boxcol="#E0FFFF",
        showP=F)
dev.off()

#输出列线图的风险打分文件
nomoRisk=predict(res.cox,data=dat1,type="risk")#进行列线图的结果预测打分
rt=cbind(dat1,Nomogram=nomoRisk)

write.csv(rt,file = "./nomoRisk.csv",row.names = F)

#校准曲线，cox回归的模型构建  注意调整B值
rt$time <- rt$time*365
pdf(file = "calibration.pdf",width =4.5 ,height = 6)
#1年校准曲线
f <- cph(Surv(time,status) ~ Nomogram,x=T,y=T,surv=T,data=rt,time.inc=1*365)#比例风险模型的拟合
cal <- calibrate(f, cmethod="KM",method="boot",u=1*365,m=(nrow(rt)/3),B=3000)
plot(cal,xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)",ylab="Observed OS (%)", lwd=1.5,col="green",sub=F)
#3年校准曲线
f <- cph(Surv(time,status) ~ Nomogram,x=T,y=T,surv=T,data=rt,time.inc=3*365)#比例风险模型的拟合
cal <- calibrate(f, cmethod="KM",method="boot",u=3*365,m=(nrow(rt)/3),B=3000)
plot(cal,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",lwd=1.5,col="blue",sub=F,add=T)
#5年校准曲线
f <- cph(Surv(time,status) ~ Nomogram,x=T,y=T,surv=T,data=rt,time.inc=5*365)#比例风险模型的拟合
cal <- calibrate(f, cmethod="KM",method="boot",u=5*365,m=(nrow(rt)/3),B=3000)
plot(cal,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",lwd=1.5,col="red",sub=F,add=T)
legend("bottomright",c('1-year','3-year','5-year'),
       col=c('green',"blue","red"),lwd=1.5,bty='n')
dev.off()

#DCA分析
library(rmda)

Age<- decision_curve(status~Age,data= rt,
                     family = binomial(link ='logit'),
                     thresholds= seq(0,1, by = 0.01),
                     confidence.intervals = 0.95,
                     study.design = 'case-control',
                     population.prevalence = 0.3)
Gender<- decision_curve(status~Gender,data= rt,
                        family = binomial(link ='logit'),
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals = 0.95,
                        study.design = 'case-control',
                        population.prevalence = 0.3)
Stage<- decision_curve(status~Stage,data= rt,
                       family = binomial(link ='logit'),
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals = 0.95,
                       study.design = 'case-control',
                       population.prevalence = 0.3)
MPCDI<- decision_curve(status~MPCDI,data= rt,
                       family = binomial(link ='logit'),
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals = 0.95,
                       study.design = 'case-control',
                       population.prevalence = 0.3)
Nomogram<- decision_curve(status~Nomogram,data= rt,
                          family = binomial(link ='logit'),
                          thresholds= seq(0,1, by = 0.01),
                          confidence.intervals = 0.95,
                          study.design = 'case-control',
                          population.prevalence = 0.3)
List<- list(Age,Gender,Stage,MPCDI,Nomogram)
pdf("DCA.pdf",width = 5,height = 5)
plot_decision_curve(List,
                    curve.names=c('Age','Gender','Stage','MPCDI','Nomogram'),
                    col= c('#008000','#00FFFF','#F4A460','#9370DB','#FF0000'),
                    confidence.intervals=FALSE,xlim = c(0,0.4))
dev.off()

#surv
surv <- rt[,c(1,2,3,8)]
library(patchwork)
surv <- read.csv("./nomoRisk.csv",header = T,row.names = 1)
surv <- surv[,c(1,2,7)]

surv$status <- ifelse(surv$status ==0,0,1)

res.cut=surv_cutpoint(surv,time = "time",event = "status",variables = c("Nomogram"))
res.cat=surv_categorize(res.cut)
surv<- res.cat
colnames(surv)[3] <- "group"

fit <- survfit(Surv(time,status) ~ group,data = surv)
res.cox <- coxph(Surv(time,status) ~group,data = surv)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ group,data = surv)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")
pdf("./surv.pdf",width = 5,height = 6)
ggsurvplot(fit,
           data = surv,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =2,
           xlab = "Overall survival(years)",
           legend.title="Nomogram Score",
           legend.labs=c("High","Low"),
           palette  = "jama",
           ggtheme = theme_bw(),
           legend="top",
           risk.table = T,
           risk.table.col = 'group',
           risk.table.y.text = F,
           ncensor.plot = T,
           censor.shape="|",censor.size=4
) + labs(caption = hr)
dev.off()

#ROC
rt <- read.csv("./nomoRisk.csv",header = T,row.names = 1)
risk=rt[,c("time","status","Nomogram")]

#绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$time,delta = risk$status,
               marker = risk$Nomogram,cause = 1,
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
load("~/GSE14520_clin.rda")
GEO1 <- GSE14520_clin[GSE14520_clin$Stage!=".",]
GEO1  <- GEO1[,-7]
colnames(GEO1)[7] <- "MPCDI"
colnames(GEO1)[1] <- "Sample"

GEO1$Stage <- factor(GEO1$Stage,levels = c("I","II","III"))
GEO1$MPCDI <- ifelse(GEO1$MPCDI == "low","Low","High")
GEO1$MPCDI <- factor(GEO1$MPCDI,levels = c("High","Low"))

res.cox2=coxph(Surv(time,status) ~Stage+MPCDI,data=GEO1)
nomoRisk2=predict(res.cox2,data=GEO1,type="risk")
rt=cbind(GEO1,Nomogram=nomoRisk2)
risk=rt[,c("time","status","Nomogram")]
#绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$time,delta = risk$status,
               marker = risk$Nomogram,cause = 1,
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
load("~/GSE116174_clin.rda")

GEO2 <- GSE116174_clin[GSE116174_clin$Stage!=".",]
GEO2  <- GEO2[,-7]
colnames(GEO2)[1] <- "Sample"
colnames(GEO2)[7] <- "MPCDI"

GEO2$Stage <- factor(GEO2$Stage,levels = c("I","II","III"))
GEO2$MPCDI <- ifelse(GEO2$MPCDI == "low","Low","High")
GEO2$MPCDI <- factor(GEO2$MPCDI,levels = c("High","Low"))

res.cox3=coxph(Surv(time,status) ~Stage+MPCDI,data=GEO2)
nomoRisk2=predict(res.cox3,data=GEO2,type="risk")
rt=cbind(GEO2,Nomogram=nomoRisk2)
risk=rt[,c("time","status","Nomogram")]
#绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$time,delta = risk$status,
               marker = risk$Nomogram,cause = 1,
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
load("~/ICGC_clin.rda")

icgc <- ICGC_clin[ICGC_clin$Stage!=".",]
icgc  <- icgc[,-7]
colnames(icgc)[1] <- "Sample"
colnames(icgc)[7] <- "MPCDI"

icgc$Stage <- ifelse(icgc$Stage == 1,"I",
                     ifelse(icgc$Stage == 2,"II",
                            ifelse(icgc$Stage==3,"III","IV")))
icgc$Stage <- factor(icgc$Stage,levels = c("I","II","III","IV"))
icgc$MPCDI <- ifelse(icgc$MPCDI == "low","Low","High")
icgc$MPCDI <- factor(icgc$MPCDI,levels = c("High","Low"))

res.cox4=coxph(Surv(time,status) ~Stage+MPCDI,data=icgc)
nomoRisk2=predict(res.cox4,data=icgc,type="risk")
rt=cbind(icgc,Nomogram=nomoRisk2)
risk=rt[,c("time","status","Nomogram")]
#绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$time,delta = risk$status,
               marker = risk$Nomogram,cause = 1,
               weighting = "aalen",
               times = c(1,3,5),ROC = TRUE)
pdf(file="ICGC-ROC.pdf",width = 5.5,height = 5.5)
plot(ROC_rt,time=1,col="#FA8072",title=F,lwd=2)
plot(ROC_rt,time=3,col="#FFD700",add=T,title=F,lwd=2)
title(main="ICGC")
legend("bottomright",
       c(paste0("AUC at 1 years: ",sprintf("%.3f",ROC_rt$AUC[1])),
         paste0("AUC at 3 years: ",sprintf("%.3f",ROC_rt$AUC[2]))),
       col=c("#FA8072","#FFD700"),lwd = 2,bty = "n")
dev.off()
