
fig=read.csv("./job_example/PDO_DEG.csv")
library(ggplot2)
library(readxl)
library(ggrepel)
fig[which(fig$log2FoldChange >= 1 & fig$padj < 0.05),'sig'] <- 'upregulated'
fig[which(fig$log2FoldChange <= -1 & fig$padj < 0.05),'sig'] <- 'downregulated'
fig[which(fig$padj %in% NA),'sig'] <- 'no diff'
fig$sig[which(is.na(fig$sig))] <- 'no diff'

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
  ) 
p
ggsave("./job_example/volcano+FGA LOW MSS VS. FGA high MSS.png",p,width=7,height=4)
