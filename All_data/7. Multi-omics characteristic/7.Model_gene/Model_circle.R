library(circlize)
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
gene_name <- read.csv("./gene_4.csv",header = T)
gene_name <- gene_name$gene
deg_gene <- gtf_data[gtf_data$gene_name %in% gene_name,]%>% distinct(gene_name,.keep_all = T)

data <- read.table("~/gencode.v22.annotation.gene.probeMap",sep = "\t" , header = T)
data <- data[data$gene %in% gene_name,]

cnv_files <- list.files(path = ".",
                        pattern = "tsv.gz",
                        full.names = T,
                        recursive = T)



tmp <- readr::read_tsv(cnv_files)
colnames(tmp)[1]<-"gene_id"

load("~/TCGA_RS.rda")
sample <- paste0(rownames(TCGA_RS),"-01A",sep="")
sample2 <- intersect(sample,colnames(tmp))
tmp <- tmp[,c("gene_id",sample2)]

tmp1 <- tmp[tmp$gene_id %in% deg_gene$gene_id,]
tmp2 <- merge(tmp1,deg_gene,by="gene_id")
tmp3 <- tmp2[,-1]
tmp4 <- tmp3 %>% dplyr::select(gene_name,everything())
rownames(tmp4) <- tmp4$gene_name
tmp5 <- tmp4[,-1]
tmp6 <- transmute(tmp5,Gain = apply(tmp5, 1, function(x){
  length(which(x == 1))
}),Loss = apply(tmp5, 1, function(x){length(which(x == -1))}))

tmp6$Gene <- rownames(tmp6)
tmp7 <- data[,c(2:5)]
colnames(tmp7)[1] <-"Gene"
tmp8 <- merge(tmp7,tmp6,by="Gene")
tmp9 <- tmp8[mixedorder(tmp8$chrom),]
tmp10 <- tmp9
tmp10$chrom <- factor(tmp10$chrom,levels = mixedsort(unique(tmp10$chrom)))
tmp10$Loss <- -tmp10$Loss
tmp11 <- tmp10[,c(2,3,4,1,5,6)]
colnames(tmp11)[1] <- "Chromosome"
write.csv(tmp11,file = "cnv_gene.csv")
tmp11 <- read.csv("./cnv_gene.csv",row.names = 1,header = T)


#SNV
load("~/laml.rda")
snv <- laml@data
snv <- snv[snv$Hugo_Symbol%in%deg_gene$gene_name,]

snv2 <- snv[,c(5,6,7,9,10,1)] %>%as.data.frame()
write.csv(snv2,file = "./snv_gene.csv")
snv3 <- snv2[,1:4]
snv3$Variant_Classification <- ifelse(snv3$Variant_Classification=="Missense_Mutation",1,2)
snv4 <- snv2[,c(1,2,3,5)]
snv4$Variant_Type <- ifelse(snv4$Variant_Type=="SNP",3,4)


pdf("./circle.pdf",width = 12,height = 12)
circos.genomicInitialize(tmp11,plotType = c("axis", "labels"),
                         track.height = 0.1)
chromosome_colors <- c("chr6"="#A13B46","chr1"="#B0CFE4","chr14"="#BEBDDF","chr22"="#F6B190")
circos.track(
  ylim = c(0, 1), 
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    col = chromosome_colors[chr]
    circos.rect(xlim[1], 0, xlim[2], 1, col = col)
    circos.text(
      mean(xlim), mean(ylim), chr, cex = 0.7,
      col = "white", facing = "inside",
      niceFacing = TRUE
    )
  }, 
  track.height = 0.1, bg.border = NA
)

circos.genomicTrack(tmp11, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 2, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "#D87070", "#95BAA6"), ...)
                      circos.genomicRect(region, value, ytop.column = 3, ybottom = 0, 
                                         col = ifelse(value[[1]] > 0, "#95BAA6", "#D87070"), ...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    },
                    track.height = 0.1
                    )
#
color_assign <- colorRamp2(breaks = c(1,2,3,4), color = c( '#84BA42', '#7ABBDB',"#682487","#DBB428"))

circos.genomicTrackPlotRegion(
  snv3, track.height = 0.1, bg.border = 'black', bg.lwd = 0.4,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, pch = 16, cex = 1.5, col = color_assign(value[[1]]), ...)
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
  } )

circos.genomicTrackPlotRegion(
  snv4, track.height = 0.1, bg.border = 'black', bg.lwd = 0.4,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, pch = 17, cex = 1.5, col = color_assign(value[[1]]), ...)
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
  } )




circos.genomicLabels(tmp11, labels.column = 4,connection_height = mm_h(2),
                     labels_height =cm_h(1.0), side = "inside")

dev.off()
