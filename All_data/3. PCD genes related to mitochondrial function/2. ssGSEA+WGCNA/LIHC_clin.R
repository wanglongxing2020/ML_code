library(tidyverse)
#TCGA-LIHC临床样本筛选
clin = read_tsv("./TCGA-LIHC_survival.tsv")
clin1 <-clin[,c("case_submitter_id","age_at_index","gender","vital_status","days_to_last_follow_up",
                "ajcc_pathologic_m","ajcc_pathologic_n",
                "ajcc_pathologic_stage","ajcc_pathologic_t")]
colnames(clin1) <- c("Tumor_Sample_Barcode","Age","Gender","status","time","M","N","Stage","T")

clin1 <- clin1[!duplicated(clin1$Tumor_Sample_Barcode),]
clin1 <- clin1[clin1$status != "'--",]
clin1$M <- gsub(pattern = "[a-d]","",clin1$M)
clin1$N <- gsub(pattern = "[a-d]","",clin1$N)
clin1$Stage <- gsub(pattern = "[A-D]","",clin1$Stage)
clin1$T <- gsub(pattern = "[a-d]","",clin1$T)


clin_ucsc <- read_tsv("./TCGA-LIHC_UCSC_survival.tsv")

clin_ucsc$group <- ifelse(as.numeric(substr(clin_ucsc$sample,14,15)) > 10,"normal","tumor")
clin_ucsc <- clin_ucsc[clin_ucsc$group == "tumor",]

clin_ucsc <- clin_ucsc[,c(3,2,4)]
colnames(clin_ucsc) <- c("Sample","status","time")

clin_tcga <- clin1[,c(1,2,3,6,7,8,9)]
colnames(clin_tcga)[1] <- "Sample"
clin_all <- merge(clin_ucsc,clin_tcga,by="Sample")
write.csv(clin_all,file = "./clin_all.csv")
