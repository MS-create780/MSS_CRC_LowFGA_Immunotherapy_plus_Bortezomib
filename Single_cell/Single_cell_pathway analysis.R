###epithelial_cell_pathway analysis###
expr_T=readRDS("./job_example/figures/epi_cluster_deleted others.rds")
library(GSVA)
library(GSEABase)
library(msigdbr)
colnames(expr_T@meta.data)
expr_T=subset(expr_T,subset = SpecimenType == "T")
result=read.csv("./job_example/epi_tumor_curated_pathway_single.csv")
rownames(result)=result[,1]
result=result[,-1]
library(limma)
meta_data=as.data.frame(expr_T@meta.data)
meta_data$ID=sub("-(.*)$", ".\\1", rownames(meta_data))
list=colnames(result)
result_metadata <- meta_data %>% filter(ID %in% list)
meta_ordered <- result_metadata %>%
  dplyr::slice(match(list, result_metadata$ID))
group <- meta_ordered[,c(35,34)]
group$subtypes[group$subtypes=="MSI-MSS"]="MSI_MSS" ##"MSI_MSS" refers to "FGA LOW MSS"
group$subtypes =as.factor(group$subtypes)
rownames(group)=group$ID
group=group[,-1,drop=F]
group=group$subtypes
desigN <- model.matrix(~ 0 + group) 
colnames(desigN) <- levels(group)
fit = lmFit(result, desigN)
cont.matrix <- makeContrasts(contrasts = c("MSI_MSS - MSS"), levels = desigN)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2,adjust='fdr', coef=1, number=Inf)
cluster2_diff <- na.omit(diff)
cluster2_diff$pathway <- rownames(cluster2_diff)
write.csv(cluster2_diff,"./job_example/epi_curated_MSI-MSS+MSS.csv")


