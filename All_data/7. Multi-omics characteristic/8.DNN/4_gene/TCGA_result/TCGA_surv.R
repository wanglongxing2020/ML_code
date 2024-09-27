library(tidyverse)
library("survival")
library("survminer")
library(patchwork)
#Entire
pred <- read.csv("./Entirepreds.csv",header = T,row.names = 1)
colnames(pred) <- "RS_pred"
pd <- read.csv("./pd_359.csv",header = T,row.names = 1)

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

pdf(file = "./TCGA_surv.pdf",width = 4,height = 6)
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


