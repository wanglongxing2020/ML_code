library(tidyverse)
load("~/TCGA_tumor_tpm.rda")
#TCGA-LIHC tpm表达矩阵
TCGA[1:4,1:4]

#临床信息
load("~/TCGA_RS.rda")
clin_info <- TCGA_RS
clin_info$RS_group <- ifelse(clin_info$RS_group=="high","High","Low")
dt <- as.data.frame(t(TCGA))
dt <- dt[rownames(clin_info),]
identical(rownames(dt), rownames(clin_info)) # 保险检查一下行名是否一致
dt$score <- clin_info$RS # 在表达矩阵中新增风险评分列
dt$riskgroup <- clin_info$RS_group # 在表达矩阵中新增风险分组列
df <- dt[order(dt$score, decreasing = F),] # 按风险评分把矩阵升序排列

df$id <- c(1:length(df$score)) # 新增id列  
df$id2 <- paste(df$riskgroup, df$id, sep = '_') # 将风险分组和id串联  
rownames(df) <- df$id2 # 修改为行名

df <- df[,1:16921] # 去掉新增列，仅保留表达矩阵  
df <- t(df) # 转置  
# 将矩阵重新转换为数值型：  
df2 <- apply(df, 2, as.numeric)  
row.names(df2) <- row.names(df)  
df2[1:6,1:6]

write.table(df2, file = 'TIDE_LIHC.txt', sep = "\t", quote = F, row.names = T) # 矩阵保存到本地