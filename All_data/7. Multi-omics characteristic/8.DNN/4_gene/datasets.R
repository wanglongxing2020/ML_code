library(tidyverse)
library(pROC)
library(ggplot2)
set.seed(123)
##训练集测试集划分1：1 每类数据集都包含一半的low和一半的high
data <- read.csv(file = "./pd_359.csv",header = T,row.names = 1)
low <- data[data$RS=="low",]
high <-data[data$RS=="high",]
train_sample1 <-sample(rownames(low),round(258/2))
test_sample1 <- rownames(low)[!rownames(low)%in%train_sample1]
length(intersect(train_sample1,test_sample1))
dat1 <- data[train_sample1,]
dat2 <- data[test_sample1,]
train_sample2 <-sample(rownames(high),51)
test_sample2 <- rownames(high)[!rownames(high)%in%train_sample2]
length(intersect(train_sample2,test_sample2))
dat1_1 <- data[train_sample2,]
dat2_2 <- data[test_sample2,]
train <- rbind(dat1,dat1_1)
test <- rbind(dat2,dat2_2)
write.csv(train,file = "./train.csv")
write.csv(test,file = "./test.csv")


#----train
#mRNA
load("~/DNN/4_gene/mRNAdata.rda")
data <- mRNA
train_sample <- read.csv("./train.csv",header = T,row.names = 1)  #注意样本量
data <- data[rownames(train_sample),]

for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
train_data <- data
write.csv(train_data,file = "./train_data.csv")

pd <- read.csv("./pd_359.csv",header = T,row.names = 1)
train_pd <- pd[rownames(train_sample),]
train_pd <- train_pd[,3,drop=F]

train_pd$RS <- ifelse(train_pd$RS=="low",0,1)
write.csv(train_pd,file = "./train_pd.csv")


#lncRNA
load("~/DNN/4_gene/lncRNAdata.rda")
data <- lncRNA
train_sample <- read.csv("./train.csv",header = T,row.names = 1)
data <- data[rownames(train_sample),]
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
train_data <- data
write.csv(train_data,file = "./train_lncdata.csv")

#Methylation
load("~/DNN/4_gene/Methydata.rda")
data <- Methy
train_sample <- read.csv("./train.csv",header = T,row.names = 1)
data <- data[rownames(train_sample),]
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}

train_data <- data
write.csv(train_data,file = "./train_methydata.csv")

#miRNA
load("~/DNN/4_gene/miRNAdata.rda")
data <- miRNA
train_sample <- read.csv("./train.csv",header = T,row.names = 1)
data <- data[rownames(train_sample),]
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}

train_data <- data
write.csv(train_data,file = "./train_midata.csv")

#-----test
#mRNA
load("~/DNN/4_gene/mRNAdata.rda")
data <- mRNA
test_sample <- read.csv("./test.csv",header = T,row.names = 1)
data <- data[rownames(test_sample),]
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}

test_data <- data
write.csv(test_data,file = "./test_data.csv")

pd <- read.csv("./pd_359.csv",header = T,row.names = 1)
test_pd <- pd[rownames(test_sample),]
test_pd <- test_pd[,3,drop=F]

test_pd$RS <- ifelse(test_pd$RS=="low",0,1)
write.csv(test_pd,file = "./test_pd.csv")



#lncRNA
load("~/DNN/4_gene/lncRNAdata.rda")
data <- lncRNA
test_sample <- read.csv("./test.csv",header = T,row.names = 1)
data <- data[rownames(test_sample),]
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
test_data <- data
write.csv(test_data,file = "./test_lncdata.csv")

#Methylation
load("~/DNN/4_gene/Methydata.rda")
data <- Methy
test_sample <- read.csv("./test.csv",header = T,row.names = 1)
data <- data[rownames(test_sample),]
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}

test_data <- data
write.csv(test_data,file = "./test_methydata.csv")

#miRNA
load("~/DNN/4_gene/miRNAdata.rda")
data <- miRNA
test_sample <- read.csv("./test.csv",header = T,row.names = 1)
data <- data[rownames(test_sample),]
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}

test_data <- data
write.csv(test_data,file = "./test_midata.csv")

#-----Entire
#mRNA
load("~/DNN/4_gene/mRNAdata.rda")
data <- mRNA
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
Entire_sample <- read.csv("./pd_359.csv",header = T,row.names = 1)
Entire_data <- data[rownames(Entire_sample),]
write.csv(Entire_data,file = "./Entire_data.csv")

pd <- read.csv("./pd_359.csv",header = T,row.names = 1)
Entire_pd <- pd[rownames(Entire_sample),]
Entire_pd <- Entire_pd[,3,drop=F]

Entire_pd$RS <- ifelse(Entire_pd$RS=="low",0,1)
write.csv(Entire_pd,file = "./Entire_pd.csv")



#lncRNA
load("~/DNN/4_gene/lncRNAdata.rda")
data <- lncRNA
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
Entire_sample <- read.csv("./pd_359.csv",header = T,row.names = 1)
Entire_data <- data[rownames(Entire_sample),]
write.csv(Entire_data,file = "./Entire_lncdata.csv")

#Methylation
load("~/DNN/4_gene/Methydata.rda")
data <- Methy
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}
Entire_sample <- read.csv("./pd_359.csv",header = T,row.names = 1)
Entire_data <- data[rownames(Entire_sample),]
write.csv(Entire_data,file = "./Entire_methydata.csv")

#miRNA
load("~/DNN/4_gene/miRNAdata.rda")
data <- miRNA
for (i in 1:100){
  data[,i] <- (data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
}

Entire_sample <- read.csv("./pd_359.csv",header = T,row.names = 1)
Entire_data <- data[rownames(Entire_sample),]
write.csv(Entire_data,file = "./Entire_midata.csv")

