library(tidyverse)
library(dplyr)
library("rtracklayer")
library("gtools")
library(reshape2)
library(stringr)
gtf_data = import('~/gencode.v22.annotation.gtf')
#提取gene_id,gene_name
gtf_data = as.data.frame(gtf_data) %>% dplyr::select(10,13)
#High-OPCDI
gene_name <- read.csv("./High_alldeg_mRNA.csv",header = T)
gene_name <- gene_name[gene_name$adj.P.Val < 0.05 & abs(gene_name$logFC) >1,]
gene_name <- gene_name$X
deg_gene <- gtf_data[gtf_data$gene_name %in% gene_name,]%>% distinct(gene_name,.keep_all = T)

data <- read.table("~/gencode.v22.annotation.gene.probeMap",sep = "\t" , header = T)
data <- data[data$gene %in% gene_name,]

load("~/TCGA_RS.rda")
TCGA_RS_group <- TCGA_RS[,-4]
colnames(TCGA_RS_group)[3] <- "RS"
High_sample <- rownames(TCGA_RS_group)[TCGA_RS_group$RS =="high"]

cnv_files <- list.files(path = ".",
                        pattern = "tsv.gz",
                        full.names = T,
                        recursive = T)



tmp <- readr::read_tsv(cnv_files)
colnames(tmp)[1]<-"gene_id"
tmp1 <- tmp[tmp$gene_id %in% deg_gene$gene_id,]
tmp2 <- merge(tmp1,deg_gene,by="gene_id")
tmp3 <- tmp2[,-1]
tmp4 <- tmp3 %>% select(gene_name,everything())
rownames(tmp4) <- tmp4$gene_name
tmp5 <- tmp4[,-1]
tmp5 <- tmp5[,substr(colnames(tmp5),15,15)=="1"]
colnames(tmp5) <- substr(colnames(tmp5),1,12)
tmp5 <- tmp5[,intersect(High_sample,colnames(tmp5))]
tmp6 <- transmute(tmp5,Gain = apply(tmp5, 1, function(x){
  length(which(x == 1))
}),Loss = apply(tmp5, 1, function(x){length(which(x == -1))}))

tmp6$Gene <- rownames(tmp6)
tmp7 <- data[,c(2,3)]
colnames(tmp7)[1] <- "Gene"
tmp8 <- merge(tmp7,tmp6,by="Gene")
tmp8 <-  tmp8 %>% distinct(Gene,.keep_all = T)
tmp9 <- tmp8[mixedorder(tmp8$chrom),]
tmp10 <- tmp9
tmp10$chrom <- factor(tmp10$chrom,levels = mixedsort(unique(tmp10$chrom)))
write.csv(tmp10,file = "High_cnv_deg.csv")

tmp11 <- tmp10
tmp11$Loss <- paste0("-",tmp11$Loss)%>%as.integer()

tmp12 <- melt(tmp11,
              id.vars =c("Gene","chrom"),
              measure.vars= c("Gain","Loss"),
              variable.name='CNV',
              value.name="Counts")



tmp12$group <- str_split_fixed(tmp12$chrom,"chr",2)[,2]
tmp12$group <-factor(tmp12$group,levels = unique(tmp12$group))
pdf("./High_cnv_deg.pdf",width =16,height = 6 )
ggplot(tmp12,aes(Gene,Counts,fill=CNV))+geom_bar(stat = "identity",position = 'dodge')+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.direction = "horizontal",
        legend.position = "top")+
  geom_hline(yintercept = 0)+
  facet_grid(cols=vars(`group`),scale='free_x',space='free_x')+
  theme(strip.background.x  = element_rect(fill = c("lightblue")))
dev.off()


#Low-OPCDI
gene_name <- read.csv("./Low_alldeg_mRNA.csv",header = T)
gene_name <- gene_name[gene_name$adj.P.Val < 0.05 & abs(gene_name$logFC) >1,]
gene_name <- gene_name$X
deg_gene <- gtf_data[gtf_data$gene_name %in% gene_name,]%>% distinct(gene_name,.keep_all = T)

data <- read.table("E~/gencode.v22.annotation.gene.probeMap",sep = "\t" , header = T)
data <- data[data$gene %in% gene_name,]

load("~/TCGA_RS.rda")
TCGA_RS_group <- TCGA_RS[,-4]
colnames(TCGA_RS_group)[3] <- "RS"
Low_sample <- rownames(TCGA_RS_group)[TCGA_RS_group$RS =="low"]

cnv_files <- list.files(path = ".",
                        pattern = "tsv.gz",
                        full.names = T,
                        recursive = T)



tmp <- readr::read_tsv(cnv_files)
colnames(tmp)[1]<-"gene_id"
tmp1 <- tmp[tmp$gene_id %in% deg_gene$gene_id,]
tmp2 <- merge(tmp1,deg_gene,by="gene_id")
tmp3 <- tmp2[,-1]
tmp4 <- tmp3 %>% select(gene_name,everything())
rownames(tmp4) <- tmp4$gene_name
tmp5 <- tmp4[,-1]
tmp5 <- tmp5[,substr(colnames(tmp5),15,15)=="1"]
colnames(tmp5) <- substr(colnames(tmp5),1,12)
sample <- intersect(colnames(tmp5),Low_sample)
tmp5 <- tmp5[,sample]
tmp6 <- transmute(tmp5,Gain = apply(tmp5, 1, function(x){
  length(which(x == 1))
}),Loss = apply(tmp5, 1, function(x){length(which(x == -1))}))

tmp6$Gene <- rownames(tmp6)
tmp7 <- data[,c(2,3)]
colnames(tmp7)[1] <- "Gene"
tmp8 <- merge(tmp7,tmp6,by="Gene")
tmp8 <-  tmp8 %>% distinct(Gene,.keep_all = T)
tmp9 <- tmp8[mixedorder(tmp8$chrom),]
tmp10 <- tmp9
tmp10$chrom <- factor(tmp10$chrom,levels = mixedsort(unique(tmp10$chrom)))
write.csv(tmp10,file = "Low_cnv_deg.csv")

tmp11 <- tmp10
tmp11$Loss <- paste0("-",tmp11$Loss)%>%as.integer()

tmp12 <- melt(tmp11,
              id.vars =c("Gene","chrom"),
              measure.vars= c("Gain","Loss"),
              variable.name='CNV',
              value.name="Counts")



tmp12$group <- str_split_fixed(tmp12$chrom,"chr",2)[,2]
tmp12$group <-factor(tmp12$group,levels = unique(tmp12$group))
pdf("./Low_cnv_deg.pdf",width =16,height = 6 )
ggplot(tmp12,aes(Gene,Counts,fill=CNV))+geom_bar(stat = "identity",position = 'dodge')+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.direction = "horizontal",
        legend.position = "top")+
  geom_hline(yintercept = 0)+
  facet_grid(cols=vars(`group`),scale='free_x',space='free_x')+
  theme(strip.background.x  = element_rect(fill = c("lightblue")))
dev.off()
