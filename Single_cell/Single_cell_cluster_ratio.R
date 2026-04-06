#epithelial cell_cluster——MSI+MSS####
Mdata=read.csv("./job_example/Epithelial_cells_number.csv",sep=",",header=T)
colnames(Mdata)
table(Mdata$ProcessingMethod)
MCD45=subset(Mdata,ProcessingMethod=="unosrted")
MCD45=as.data.frame(MCD45)
MCD451=aggregate(x=MCD45$n,by=list(MCD45$PatientBarcode,MCD45$SpecimenType,MCD45$newsite,MCD45$MSI.state),FUN=sum)
colnames(MCD451)[1]="PatientBarcode"
colnames(MCD451)[2]="SpecimenType"
colnames(MCD451)[3]="newsite"
colnames(MCD451)[4]="MSI.state"
colnames(MCD45)[1]="a"
CD4ratio=merge(MCD45,MCD451,how='outter',on='PatientBarcode',all=T)
CD4ratio=CD4ratio[,-5]
CD4ratio$ratio=CD4ratio$n/CD4ratio$x
CD4ratio$SpecimenType[CD4ratio$SpecimenType!="N"]="T"
CD4ratio$MSI.state[CD4ratio$MSI.state!="MSS"]="MSI"
#CD4ratio=subset(CD4ratio,CD4ratio$MSI.state!="MSS")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
CD4ratio$SpecimenType=as.factor(CD4ratio$SpecimenType)
CD4ratio$seurat_clusters=as.factor(CD4ratio$seurat_clusters)
a_list=c(11,13,14,18,19,20)
CD4ratio=CD4ratio %>% group_by(seurat_clusters) %>%
  filter(!seurat_clusters %in% a_list) %>%
  ungroup()
p <- ggboxplot(CD4ratio, x = "seurat_clusters", y = "ratio",size=0.5, bxp.errorbar = T, 
               color = "MSI.state",  palette =c( "#2dabb2", "#daab36"),
               add.params = list(size=0.5), outlier.shape = NA, 
               add = "point") + 
  stat_compare_means(aes(group = MSI.state), label = "p.signif",hide.ns = TRUE)
p
ggsave("./job_example/MSS+MSI-epi-ratio_compare.pdf",p,height=6,width=10)


##epithelial cluster——N+T########
Mdata=read.csv("./job_example/Epithelial_cells_number.csv",sep=",",header=T)
colnames(Mdata)
table(Mdata$ProcessingMethod)
MCD45=subset(Mdata,ProcessingMethod=="unosrted")
MCD45=as.data.frame(MCD45)
MCD45$SpecimenType[MCD45$SpecimenType=="TA"]="T"
MCD45$SpecimenType[MCD45$SpecimenType=="TB"]="T"
MCD451=aggregate(x=MCD45$n,by=list(MCD45$PatientBarcode,MCD45$SpecimenType,MCD45$newsite,MCD45$MSI.state),FUN=sum)
colnames(MCD451)[1]="PatientBarcode"
colnames(MCD451)[2]="SpecimenType"
colnames(MCD451)[3]="newsite"
colnames(MCD451)[4]="MSI.state"
colnames(MCD45)[1]="a"
CD4ratio=merge(MCD45,MCD451,how='outter',on='PatientBarcode',all=T)
CD4ratio=CD4ratio[,-5]
CD4ratio$ratio=CD4ratio$n/CD4ratio$x
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
CD4ratio$newsite=as.factor(CD4ratio$newsite)
CD4ratio$seurat_clusters=as.factor(CD4ratio$seurat_clusters)
a_list=c(11,13,14,18,19,20)
CD4ratio=CD4ratio %>% group_by(seurat_clusters) %>%
  filter(!seurat_clusters %in% a_list) %>%
  ungroup()
p <- ggboxplot(CD4ratio, x = "seurat_clusters", y = "ratio",size=0.5, bxp.errorbar = T, 
               color = "SpecimenType",  palette =c( "#2dabb2", "#daab36"),
               add.params = list(size=0.5), outlier.shape = NA, 
               add = "point") + 
  stat_compare_means(aes(group = SpecimenType), label = "p.signif",hide.ns = TRUE)
p
ggsave("./job_example/T+N-epi-ratio_compare.pdf",p,height=6,width=10)