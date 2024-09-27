
load("~/DNN/4_gene/mRNA/DNN_100mRNA_expr.rda")
load("~/DNN/4_gene/miRNA/DNN_100miRNA_expr.rda")
load("~/DNN/4_gene/lncRNA/DNN_100lncRNA_expr.rda")
load("~/DNN/4_gene/Methylation/DNN_100Methylation_expr.rda")

Methy <- t(DNN_100Methylation_expr)
mRNA <- DNN_100mRNA_expr
miRNA <- DNN_100miRNA_expr
lncRNA <- DNN_100lncRNA_expr

rownames(Methy) <-substr(rownames(Methy),1,12)
rownames(miRNA) <-substr(rownames(miRNA),1,12)
rownames(lncRNA) <-substr(rownames(lncRNA),1,12)
rownames(mRNA) <- substr(rownames(mRNA),1,12)
#359个样本
Sample <- Reduce(intersect, list(rownames(Methy), rownames(miRNA), 
                                 rownames(lncRNA), rownames(mRNA)))
Methy <- Methy[Sample,]
mRNA <- mRNA[Sample,]
lncRNA <- lncRNA[Sample,]
miRNA <- miRNA[Sample,]

save(mRNA,file = "mRNAdata.rda")
save(miRNA,file = "miRNAdata.rda")
save(lncRNA,file = "lncRNAdata.rda")
save(Methy,file = "Methydata.rda")
load("~/TCGA_RS.rda")

pd <- TCGA_RS[Sample,]
pd <- pd[,-4]
colnames(pd)[3] <- "RS"
write.csv(pd,file = "pd_359.csv")
