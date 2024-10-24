## This is a prototype of differential expression analysis results
## from DESeq2 that I used for making volcano plots for all 
## differentially expressed genes as well as genes closest to a ligand binding site 
## from ChIP sequencing for potential receptor regulated genes


e2=read.table(file="DMSO_vs_E2_DESEq2_Diff_exp_analysis.txt",header = T,sep="\t",stringsAsFactors = F)
adjp=0.01

e2$baseMean_log=log2(e2$baseMean+1)


library(dplyr)
library(biomaRt)


library(org.Hs.eg.db)

e2$Symbol <- mapIds(org.Hs.eg.db, e2$Ensembl,keytype="ENSEMBL", column="SYMBOL")
e2$Neg_log_p_val=-log10(e2$padj)



#all_combined=rbind(up,down,no_change)
library(ggrepel)

e2_sig_up=subset(e2,e2$log2FoldChange >1 & e2$padj < 0.01)
e2_sig_dn=subset(e2,e2$log2FoldChange < -1 & e2$padj < 0.01)

e2_sig=rbind(e2_sig_up,e2_sig_dn)

e2_sig_up$e2_Direction=ifelse(e2_sig_up$log2FoldChange > 0,"Up")
e2_sig_dn$e2_Direction=ifelse(e2_sig_dn$log2FoldChange < 0,"Down")

e2_sig_up=e2_sig_up[order(-e2_sig_up$Neg_log_p_val),]
e2_top_10_up=e2_sig_up[1:10,]

e2_sig_dn=e2_sig_dn[order(-e2_sig_dn$Neg_log_p_val),]
e2_top_10_dn=e2_sig_dn[1:10,]


library(ggplot2)
p=ggplot(e2, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=e2, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = e2_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = e2_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")

p3 = p2+ggtitle("Volcano plot for DMSO vs E2")+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate(geom="text", x=12, y=-3, label="518 genes up after E2",color="red",size=6)+
  annotate(geom="text", x=12, y=-4, label="348 genes down after E2",color="blue",size=6)+
  xlab("Log2FoldChange")+xlim(-8,8)+ylab("-log10(p value)")+ylim(0,300)+
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf,data=e2_top_10_up,aes(label=Symbol, size = 12, color = NULL))+
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf,data=e2_top_10_dn,aes(label=Symbol, size = 12, color = NULL))+
  theme(axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 12),axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 12))
p3

tiff("DMSO vs E2 Volcano plot.tiff",width=800)
p3
dev.off()




e2_er=read.table(file="E2_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
colnames(e2_er)[1:3]=c("seqnames","start","end")
chr_keep=1:22
chr_all=c(chr_keep,"X")

e2_er=e2_er[e2_er$seqnames %in% chr_all,]

library(ChIPseeker)
library(ChIPpeakAnno)
library(AnnotationDbi)

e2_er=toGRanges(e2_er[,c(1:3)],format="BED",header=TRUE)

txdb=loadDb("hg38_transcription_txdb.Rdata")


peakAnno.e2.er <- annotatePeak(e2_er, tssRegion=c(-1000, 1000),
                               TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.e2.er=as.data.frame(peakAnno.e2.er)


e2_up_dn_er=intersect(e2_sig$Ensembl,peakAnno.e2.er$geneId)

`%ni%`=Negate(`%in%`)
e2_up_dn_not_er_target=e2_sig$Ensembl[e2_sig$Ensembl %ni% peakAnno.e2.er$geneId]

e2_up_dn_er_tab=subset(e2_sig,e2_sig$Ensembl %in% e2_up_dn_er)
e2_up_dn_not_er_tab=subset(e2_sig,e2_sig$Ensembl %in% e2_up_dn_not_er_target)


e2_up_er_tab=subset(e2_up_dn_er_tab,e2_up_dn_er_tab$log2FoldChange > 0)
e2_dn_er_tab=subset(e2_up_dn_er_tab,e2_up_dn_er_tab$log2FoldChange < 0)

e2_up_not_er_tab=subset(e2_up_dn_not_er_tab,e2_up_dn_not_er_tab$log2FoldChange > 0)
e2_dn_not_er_tab=subset(e2_up_dn_not_er_tab,e2_up_dn_not_er_tab$log2FoldChange < 0)




p=ggplot(e2, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=e2, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = e2_up_er_tab, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = e2_dn_er_tab, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")
p3 <- p2 +  geom_point(data = e2_up_not_er_tab, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="purple")
p4 <- p3 +  geom_point(data = e2_dn_not_er_tab, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="green")

p5 = p4+ggtitle("Volcano plot for DMSO vs E2")+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate(geom="text", x=-5, y=280, label="339 genes up closest to ERa",color="red",size=6)+
  annotate(geom="text", x=-5, y=270, label="193 genes down closest to ERa",color="blue",size=6)+
  annotate(geom="text", x=-5, y=260, label="179 genes up not closest to ERa",color="purple",size=6)+
  annotate(geom="text", x=-5, y=250, label="155 genes down not closest to ERa",color="green",size=6)+
  xlab("Log2FoldChange")+xlim(-8,8)+ylab("-log10(p value)")+ylim(0,300)+
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf,data=e2_top_10_up,aes(label=Symbol, size = 12, color = NULL))+
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf,data=e2_top_10_dn,aes(label=Symbol, size = 12, color = NULL))+
  theme(axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 12),axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 12))
p5






