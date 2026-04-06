
fig=read.csv("./job_example/PDO_DEG.csv")
library(ggplot2)
library(readxl)
library(ggrepel)
fig[which(fig$log2FoldChange >= 1 & fig$padj < 0.05),'sig'] <- 'upregulated'
fig[which(fig$log2FoldChange <= -1 & fig$padj < 0.05),'sig'] <- 'downregulated'
fig[which(fig$padj %in% NA),'sig'] <- 'no diff'
fig$sig[which(is.na(fig$sig))] <- 'no diff'

up_top10 <- fig[fig$sig == "upregulated", ]
up_top10 <- up_top10[order(-up_top10$log2FoldChange), ][1:10, ]


down_top10 <- fig[fig$sig == "downregulated", ]
down_top10 <- down_top10[order(down_top10$log2FoldChange), ][1:10, ]


label_genes <- rbind(up_top10, down_top10)
qq=label_genes$X
qq=qq[-c(6,18)]
qq
colnames(fig)
cut_off_qvalue = 0.05
cut_off_logFC = 1
p=ggplot(fig, aes(x = log2FoldChange, y = -log10(padj), colour=sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + xlim(c(-16, 16)) + 
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_qvalue),
             lty=4,col="black",lwd=0.8) +
  labs(x="Fold Change", y="-log10 (padj)") +
  theme_bw() +
  ggtitle("MSS with low FGA vs. Classic MSS") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank() 
  ) +  
  geom_text_repel(
    data = fig[fig$X %in% qq,],
    aes(label = X), 
    size = 3,
    point.padding = unit(1.5, "lines"),  
    box.padding = unit(1, "lines"),     
    max.overlaps = 20,  
    force = 2,           
    segment.color = "gray50",  
    segment.size = 0.3,        
    bg.color = "white",        
    bg.r = 0.1,                
    show.legend = FALSE
  )
p
ggsave("./job_example/volcano+FGA LOW MSS VS. CLASSIC MSS.png",p,width=7,height=4)
