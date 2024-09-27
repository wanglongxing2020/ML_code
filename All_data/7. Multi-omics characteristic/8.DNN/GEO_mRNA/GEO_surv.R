library(tidyverse)
library("survival")
library("survminer")
#GSE14520
pred <- read.csv("./GSE14520preds.csv",header = T,row.names = 1)
colnames(pred) <- "RS_pred"

load("~/GSE14520_RS.rda")
GSE14520_RS_group <- GSE14520_RS[,-4]
colnames(GSE14520_RS_group)[3] <- "RS"
pd <- GSE14520_RS_group

data <- cbind(pd,RS_pred=pred$RS_pred)

data$RS_pred <- ifelse(data$RS_pred >0.5,"High","Low")
fit <- survfit(Surv(time, status) ~ RS_pred, data = data)
print(fit)
res.cox <- coxph(Surv(time,status) ~RS_pred,data = data)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS_pred,data = data)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")

pdf(file = "./GSE14520_surv.pdf",width = 4,height = 6)
ggsurvplot(fit,
           data = data,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =2,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama",
           ggtheme = theme_bw(),
           legend="top",
           risk.table = T,
           risk.table.col = 'RS_pred',
           risk.table.y.text = F
) + labs(caption = hr)
dev.off()

#GSE116174
pred <- read.csv("./GSE116174preds.csv",header = T,row.names = 1)
colnames(pred) <- "RS_pred"

load("~/GSE116174_RS.rda")
GSE116174_RS_group <- GSE116174_RS[,-4]
colnames(GSE116174_RS_group)[3] <- "RS"
pd <- GSE116174_RS_group

data <- cbind(pd,RS_pred=pred$RS_pred)

data$RS_pred <- ifelse(data$RS_pred >0.5,"High","Low")
fit <- survfit(Surv(time, status) ~ RS_pred, data = data)
print(fit)
res.cox <- coxph(Surv(time,status) ~RS_pred,data = data)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS_pred,data = data)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")

pdf(file = "./GSE116174_surv.pdf",width = 4,height = 6)
ggsurvplot(fit,
           data = data,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =2,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama",
           ggtheme = theme_bw(),
           legend="top",
           risk.table = T,
           risk.table.col = 'RS_pred',
           risk.table.y.text = F
) + labs(caption = hr)
dev.off()

#ICGC
pred <- read.csv("./ICGCpreds.csv",header = T,row.names = 1)
colnames(pred) <- "RS_pred"

load("~/ICGC_RS.rda")
ICGC_RS_group <- ICGC_RS[,-4]
colnames(ICGC_RS_group)[3] <- "RS"
pd <- ICGC_RS_group

data <- cbind(pd,RS_pred=pred$RS_pred)

data$RS_pred <- ifelse(data$RS_pred >0.5,"High","Low")
fit <- survfit(Surv(time, status) ~ RS_pred, data = data)
print(fit)
res.cox <- coxph(Surv(time,status) ~RS_pred,data = data)
summary(res.cox)
data.survdiff <- survdiff(Surv(time, status) ~ RS_pred,data = data)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
hr <- paste0("HR = ",round(HR,2),"(",round(low95,2),"-",round(up95,2),")")

pdf(file = "./ICGC_surv.pdf",width = 4,height = 6)
ggsurvplot(fit,
           data = data,
           conf.int = T,
           conf.int.style="step",
           pval = T,
           pval.method = T,
           break.x.by =2,
           xlab = "Overall survival(years)",
           legend.title="MPCDI",
           legend.labs=c("High","Low"),
           palette  = "jama",
           ggtheme = theme_bw(),
           legend="top",
           risk.table = T,
           risk.table.col = 'RS_pred',
           risk.table.y.text = F
) + labs(caption = hr)
dev.off()
