#DEG for organoids data####
setwd("./job_example/")
exp=read.csv("RNA-seq_raw.csv")
gro=read.csv("FGA_organoid_group.csv")
exp=exp[,-1]
exp=exp[rowSums(exp[,c(2,17)])!=0,]
library(limma)
exp=avereps(exp)
rownames(exp)=exp[,1]
exp=as.data.frame(exp)
exp=exp[,-1]
gro=subset(gro,gro$group!="MSI")
exp2=exp[,c(1,2)]
list=gro$ID
for(i in list){
  a=exp[,colnames(exp)==i]
  exp2=cbind(exp2,a)
}
exp2=exp2[,-c(1,2)]
colnames(exp2)=list
rownames(gro)=gro[,1]
gro=gro[,-1,drop=F]
gro2<-factor(t(gro$group))
colnames(gro)[1]="group"
gro$group<-relevel(gro$group,ref = 'FGA high')###'classic MSS' refers to high FGA MSS samples
exp2=as.data.frame(lapply(exp2,as.numeric))
rownames(exp2)=rownames(exp)

library("DESeq2") 
library("dplyr")
qualified_genes <- c()
for (genes_in_sheet in rownames(exp2)) {
  qualification <- exp2[genes_in_sheet,] <= 1
  if (sum(qualification) < 0.8*length(exp2)) {
    qualified_genes <- append(qualified_genes, genes_in_sheet)
  }
}
mRNA_expr_for_DESeq <- exp2[qualified_genes,]
count=as.matrix(mRNA_expr_for_DESeq)
any(count < 0)
count[count<0] <- 0
sum(is.na(count))
which(is.na(count))
mode(count) <- "integer"
count=count[complete.cases(count),]
dds<- DESeqDataSetFromMatrix(count, colData =gro, design = ~ group)
dds_DE<-DESeq(dds)
res_DE <- results(dds_DE,contrast=c("group","FGA low","FGA high"))
res_DE=as.data.frame(res_DE)
write.csv(res_DE, 'DEseq+FGA low+FGA high.csv')
