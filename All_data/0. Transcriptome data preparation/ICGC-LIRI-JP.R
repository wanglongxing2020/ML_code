library(tidyverse)
#临床数据
clin <- read.csv("./donor.LIRI-JP.tsv",na="NA",sep = "\t")
clin <- clin[complete.cases(clin$donor_survival_time),]
clin1 <- clin[,c(1,5,6,9,15,17)]

clin1$Age <- clin1$donor_age_at_diagnosis
clin1$Gender <- clin1$donor_sex
clin1$status <- ifelse(clin1$donor_vital_status=="alive",0,1)
clin1$time <- clin1$donor_survival_time
clin1$time <- clin1$time/365
clin1$Stage <- clin1$donor_tumour_stage_at_diagnosis
clin2 <- clin1[,c(1,7:11)]

ICGC_clin <- clin2
save(ICGC_clin,file = "./ICGC_clin.rda")



specimen_file <- read_delim(file="specimen.LIRI-JP.tsv",delim="\t",col_names = TRUE)

#C:273  N:468
sepcimen_simplify_file <-   specimen_file %>% mutate(specimen_type=ifelse(specimen_type=="Primary tumour - solid tissue","C","N")) %>%
  dplyr::select(1,5,7)

sepcimen_simplify_file <-sepcimen_simplify_file[sepcimen_simplify_file$specimen_type=="C",]


exp_file <- read_delim(file = "exp_seq.LIRI-JP.tsv",delim = "\t",col_names = TRUE)
exp_file <- exp_file[,c(1,3,8,9)]

## 将两文件合并，形成合并的表达矩阵
merge_df <- merge(exp_file,sepcimen_simplify_file,by="icgc_specimen_id")


## gather 格式的表达数据
tmp <- merge_df[,c(1,5,6,3,4)] %>% tbl_df() %>% 
  unite(sample_id,icgc_specimen_id,icgc_donor_id.y,specimen_type,sep="-") 
tmp <- tmp[complete.cases(tmp),]

## 查看捐赠者数目
length(unique(exp_file$icgc_donor_id))    # [1] 232
## 查看sample 数目
length(unique(exp_file$icgc_specimen_id)) # [1] 445

## spread 格式表达数据
## 解决spread报错：https://www.jianshu.com/p/de03346a584a

save(tmp,file = "tmp.rda")

exp_matrix <-spread(tmp,key = "sample_id",value = "normalized_read_count")

exp_matrix <- exp_matrix[complete.cases(exp_matrix),]

exp_matrix <- as.data.frame(exp_matrix)
rownames(exp_matrix) <- exp_matrix$gene_id
exp_matrix <- exp_matrix[,-1]
ICGC_exp <- exp_matrix
ICGC_exp <- log2(ICGC_exp+1)
save(ICGC_exp,file = "ICGC_exp.rda")

exp <- ICGC_exp
colnames(exp) <- sapply(strsplit(colnames(exp),"\\-"),'[',2)

exp <- t(exp)
clin <- ICGC_clin
rownames(clin) <- clin$icgc_donor_id
clin <- clin[,-1]
exp_clin <- merge(clin,exp,by="row.names")
rownames(exp_clin) <- exp_clin$Row.names
exp_clin <- exp_clin[,-1]
ICGC_exp_clin <- exp_clin
save(ICGC_exp_clin,file = "ICGC_exp_clin.rda")
