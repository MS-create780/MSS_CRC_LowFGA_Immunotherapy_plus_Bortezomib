library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(ggExtra)
library(patchwork)
library(cowplot)
library(reshape2)
library(pheatmap)
library(ggsci)
library(ggpubr)
library(ggsignif)
library(ggstatsplot)
library(harmony)
expr=readRDS("./job_example/pbmc_85.rds")
meta_data=expr@meta.data
table(expr@meta.data[["Idents"]])
expr=NormalizeData(expr,normalization.method = "LogNormalize")
expr<- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000, verbose = F)
expr <- ScaleData(expr, verbose = FALSE)
expr=RunPCA(expr, verbose = FALSE)
expr <- FindNeighbors(expr, reduction = "pca", dims = 1:10, verbose = F)
expr <- RunTSNE(object = expr, dims = 1:30, reduction = "pca")
expr <- RunUMAP(expr,dims=1:20)
saveRDS(expr,"./job_example/singlecell_afterUMAP.rds")
###marker gene###3
expr <- SetIdent(expr, value = "Idents")
all_markers <- FindAllMarkers(expr,only.pos = FALSE,            
                                  min.pct = 0.25,              
                                  logfc.threshold = 0.25)
write.csv(all_markers,"./job_example/cluster_markers.csv")
###annotation###
meta_data=expr@meta.data
meta_data$cell_type=sub("([a-zA-Z]+)\\d+","\\1",meta_data$Idents)
table(meta_data$cell_type)
meta_data$cell_type[meta_data$cell_type=="B"]="B_cells"
meta_data$cell_type[meta_data$cell_type=="En"]="Endothelial_cells"
meta_data$cell_type[meta_data$cell_type=="Ep"]="Epithelial_cells"
meta_data$cell_type[meta_data$cell_type=="F"]="Fibroblast"
meta_data$cell_type[meta_data$cell_type=="M"]="Mast_cells"
meta_data$cell_type[meta_data$cell_type=="My"]="Myeloid_cells"
meta_data$cell_type[meta_data$cell_type=="P"]="Plasma_cells"
meta_data$cell_type[meta_data$cell_type=="T"]="T_cells"
colnames(meta_data)[24]="all_cell_type"
expr=AddMetaData(expr, metadata = meta_data)

#UMAP plot####
p=DimPlot(expr, reduction = "umap", label = TRUE, pt.size = 1.2, group.by ="all_cell_type")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
p2 <- p + scale_color_npg() + labs(title = "npg", tag = "B")
p2
ggsave("./job_example/all_UMAP_plot.pdf",p2,width=9,height=6)

#marker gene plot####
cols = c("gray", "coral2")
plot9 <- FeaturePlot(expr, features = c('PTPRC',"CD19","CD79A","CD3D","CD3E",
                                        "CD14","CD68","TPSAB1","KIT","PTPRB","KDR",
                                        "AQP1","PLVAP","EPCAM","KRT8","KRT18",
                                        "WFDC2","COL1A1","COL1A2","COL6A1","COL3A1",
                                        "CD14","CD33","FCGR1A","FCGR3B","CSF3R"),cols = cols, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
ggsave("./job_example/marker_feature plot.pdf",plot9,width=20,height=20)

