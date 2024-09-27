library(tidyverse)
#GSE14520
data <-read.csv("./GSE14520_data.csv",header = T,row.names = 1)

for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
load("~/GSE14520_RS.rda")
GSE14520_RS_group <- GSE14520_RS[,-4]
colnames(GSE14520_RS_group)[3] <- "RS"
vali_sample <- GSE14520_RS_group
vali_data <- data[rownames(vali_sample),]
write.csv(vali_data,file = "./GSE14520_data.csv")

pd <- GSE14520_RS_group
vali_pd <- pd[rownames(vali_sample),]
vali_pd <- vali_pd[,3,drop=F]

vali_pd$RS <- ifelse(vali_pd$RS=="low",0,1) #low:165 high:56
write.csv(vali_pd,file = "./GSE14520_pd.csv")

#GSE116174
data <-read.csv("./GSE116174_data.csv",header = T,row.names = 1)

for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
load("~/GSE116174_RS.rda")
GSE116174_RS_group <- GSE116174_RS[,-4]
colnames(GSE116174_RS_group)[3] <- "RS"
vali_sample <- GSE116174_RS_group
vali_data <- data[rownames(vali_sample),]
write.csv(vali_data,file = "./GSE116174_data.csv")

pd <- GSE116174_RS_group
vali_pd <- pd[rownames(vali_sample),]
vali_pd <- vali_pd[,3,drop=F]

vali_pd$RS <- ifelse(vali_pd$RS=="low",0,1)  #low:51 high:13
write.csv(vali_pd,file = "./GSE116174_pd.csv")

#ICGC
data <-read.csv("./ICGC_data.csv",header = T,row.names = 1)

for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
load("~/ICGC_RS.rda")
ICGC_RS_group <- ICGC_RS[,-4]
colnames(ICGC_RS_group)[3] <- "RS"
vali_sample <- ICGC_RS_group
vali_data <- data[rownames(vali_sample),]
write.csv(vali_data,file = "./ICGC_data.csv")

pd <- ICGC_RS_group
vali_pd <- pd[rownames(vali_sample),]
vali_pd <- vali_pd[,3,drop=F]

vali_pd$RS <- ifelse(vali_pd$RS=="low",0,1)  #low:90 high:141
write.csv(vali_pd,file = "./ICGC_pd.csv")

