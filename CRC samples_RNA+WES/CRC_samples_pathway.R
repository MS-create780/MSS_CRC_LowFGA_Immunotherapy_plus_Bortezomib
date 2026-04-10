#pathway analysis####
RNA=read.csv("./job_example/FGA0.34+DEseq+low+high.csv")
diff=RNA
library(enrichplot)
library(clusterProfiler)
library(data.table)
library(stringr)
library(ggplot2)
geneList <- diff$log2FoldChange
names(geneList) = diff$gene
geneList = sort(geneList, decreasing = TRUE)
geneList = geneList[!is.infinite(geneList)]
library(clusterProfiler)
hallmarks <- read.gmt("./job_example/c2.all.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList,TERM2GENE =hallmarks,pvalueCutoff=1,eps = 0,nPermSimple=10000)
head(y)
y2=as.data.frame(y)
write.csv(y,"./job_example/FGA0.34+RNA+LOW VS HIGH+hallmarker.csv")
