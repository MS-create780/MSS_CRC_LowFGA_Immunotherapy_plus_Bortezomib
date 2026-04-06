#Differential gene_DEG####
data=read.csv("./job_example/readcount.csv")
group=read.csv("./job_example/clinical.csv",header=F,encoding='UTF-8')
colnames(group)=group[1,]
group=group[-1,]
group=subset(group,!is.na(group$I_number))
group=group[-41,]
group$group[group$FGA<=0.35]="low"
group$group[group$FGA>0.35]="high"
group=subset(group,MSI.state=="MSS")
group=group[,c(4,80)]
p=group[,1]
group=as.data.frame(group[,-1])
rownames(group)=p
colnames(group)="group"

data2=as.data.frame(data[,1])
list=rownames(group)
for(i in list){
  a=data[,colnames(data)==i]
  data2=cbind(data2,a)
}
colnames(data2)[c(2:37)]=list
colnames(data2)[1]="gene"
library(limma)
data2=avereps(data2,data2$gene)
data2=as.data.frame(data2)
rownames(data2)=data2[,1]
data2=data2[,-1]
q=rownames(data2)
data2=as.data.frame(lapply(data2,as.numeric))
rownames(data2)=q
is.numeric(data2$I18)
library("DESeq2") 
library("dplyr")
qualified_genes <- c()
for (genes_in_sheet in rownames(data2)) {
  qualification <- data2[genes_in_sheet,] <= 1
  if (sum(qualification) < 0.8*length(data2)) {
    qualified_genes <- append(qualified_genes, genes_in_sheet)
  }
}
mRNA_expr_for_DESeq <- data2[qualified_genes,]
count=as.matrix(mRNA_expr_for_DESeq)
any(count < 0)
count[count<0] <- 0
sum(is.na(count))
which(is.na(count))
mode(count) <- "integer"
count=count[complete.cases(count),]

dds<- DESeqDataSetFromMatrix(count, colData =group, design = ~ group)
dds_DE<-DESeq(dds)
res_DE <- results(dds_DE,contrast=c("group","low","high"))
write.csv(res_DE, './job_example/FGA0.35+DEseq+low+high.csv')