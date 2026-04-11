diff=read.csv("PDO_DEG.csv")
diff$X=gsub("[.]","-",diff$X)
library(enrichplot)
library(clusterProfiler)
library(data.table)
library(stringr)
library(ggplot2)
colnames(diff)[1]="gene"
geneList <- diff$log2FoldChange
names(geneList) = diff$gene
geneList = sort(geneList, decreasing = TRUE)
geneList = geneList[!is.infinite(geneList)]
library(clusterProfiler)
hallmarks <- read.gmt("C:/COAD/non transcript/c2.all.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList,TERM2GENE =hallmarks,pvalueCutoff=1,eps = 0)
head(y)
y2=as.data.frame(y)
write.csv(y,"PDO_GSEA_pathway.csv")
