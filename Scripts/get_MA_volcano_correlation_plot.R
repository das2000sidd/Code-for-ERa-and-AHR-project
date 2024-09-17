setwd("~/Desktop/PhD_Project_related/EO771_WT_AhrKOmice_tumors")


norm_counts=read.table(file="CPM_normalised_counts.txt",header = T,stringsAsFactors = F,sep="\t")
wt_E0771_vs_PY8119=read.csv(file="res_wt_E0771_vs_PY8119_EO771_WT_AhrKOmice_tumors.csv",header = T,stringsAsFactors = F,sep=",")
adjp=0.01

#wt_E0771_vs_PY8119=wt_E0771_vs_PY8119[complete.cases(wt_E0771_vs_PY8119),]
wt_E0771_vs_PY8119$baseMean_log=log2(wt_E0771_vs_PY8119$baseMean+1)

library(org.Mm.eg.db)

#wt_E0771_vs_PY8119$Entrez <- mapIds(org.Mm.eg.db, wt_E0771_vs_PY8119$Ensembl,keytype="ENSEMBL", column="ENTREZID")
#wt_E0771_vs_PY8119$Symbol <- mapIds(org.Mm.eg.db, wt_E0771_vs_PY8119$Entrez,keytype="ENTREZID", column="SYMBOL")
#wt_E0771_vs_PY8119$Genename <- mapIds(org.Mm.eg.db, wt_E0771_vs_PY8119$Entrez,keytype="ENTREZID", column="GENENAME")

#wt_E0771_vs_PY8119$Entrez=as.character(wt_E0771_vs_PY8119$Entrez)
#wt_E0771_vs_PY8119$Symbol=as.character(wt_E0771_vs_PY8119$Symbol)
#wt_E0771_vs_PY8119$Genename=as.character(wt_E0771_vs_PY8119$Genename)


wt_E0771_vs_PY8119_noGm_riken=wt_E0771_vs_PY8119
#ahrko_noGm_riken=ahrko[- grep("predicted",ahrko$Gene_Name),]
#wt_E0771_vs_PY8119_noGm_riken=wt_E0771_vs_PY8119_noGm_riken[complete.cases(wt_E0771_vs_PY8119_noGm_riken),]
#wt_E0771_vs_PY8119_noGm_riken=wt_E0771_vs_PY8119_noGm_riken[- grep("predicted",wt_E0771_vs_PY8119_noGm_riken$Genename),]
#wt_E0771_vs_PY8119_noGm_riken=wt_E0771_vs_PY8119_noGm_riken[- grep("Riken",wt_E0771_vs_PY8119_noGm_riken$Genename),]
#wt_E0771_vs_PY8119_noGm_riken=wt_E0771_vs_PY8119_noGm_riken[- grep("pseudogene",wt_E0771_vs_PY8119_noGm_riken$Genename),]



library(dplyr)
library(biomaRt)

library(ggrepel)

wt_E0771_vs_PY8119_noGm_riken$Neg_log_p_val=-log10(wt_E0771_vs_PY8119_noGm_riken$padj)

wt_E0771_vs_PY8119_sig_up=subset(wt_E0771_vs_PY8119_noGm_riken,wt_E0771_vs_PY8119_noGm_riken$log2FoldChange>1 & wt_E0771_vs_PY8119_noGm_riken$padj < 0.01)
wt_E0771_vs_PY8119_sig_dn=subset(wt_E0771_vs_PY8119_noGm_riken,wt_E0771_vs_PY8119_noGm_riken$log2FoldChange < -1 & wt_E0771_vs_PY8119_noGm_riken$padj < 0.01)

wt_E0771_vs_PY8119_sig_up$wt_E0771_vs_PY8119_Direction="wt_E0771_vs_PY8119_Up"
wt_E0771_vs_PY8119_sig_dn$wt_E0771_vs_PY8119_Direction="wt_E0771_vs_PY8119_Down"

wt_E0771_vs_PY8119_sig=rbind(wt_E0771_vs_PY8119_sig_up,wt_E0771_vs_PY8119_sig_dn)

wt_E0771_vs_PY8119_sig=wt_E0771_vs_PY8119_sig[,c("Ensembl","wt_E0771_vs_PY8119_Direction")]

#wt_E0771_vs_PY8119_noGm_riken$Ensembl=rownames(wt_E0771_vs_PY8119_noGm_riken)

wt_E0771_vs_PY8119_noGm_riken=left_join(wt_E0771_vs_PY8119_noGm_riken,wt_E0771_vs_PY8119_sig,by=c("Ensembl"))
wt_E0771_vs_PY8119_noGm_riken$wt_E0771_vs_PY8119_Direction[is.na(wt_E0771_vs_PY8119_noGm_riken$wt_E0771_vs_PY8119_Direction)]="No sig change"




#table(wt_E0771_vs_PY8119=wt_E0771_vs_PY8119_sig_up$wt_E0771_vs_PY8119_Direction,`wt_E0771_vs_PY8119 Up`=wt_E0771_vs_PY8119_sig_up$wt_E0771_vs_PY8119_Direction)


wt_E0771_vs_PY8119_noGm_riken$baseMean_log=log2(wt_E0771_vs_PY8119_noGm_riken$baseMean+1)

wt_E0771_vs_PY8119_noGm_riken_inf=wt_E0771_vs_PY8119_noGm_riken[grep("Inf",wt_E0771_vs_PY8119_noGm_riken$Neg_log_p_val),]
wt_E0771_vs_PY8119_noGm_riken_not_inf=wt_E0771_vs_PY8119_noGm_riken[- grep("Inf",wt_E0771_vs_PY8119_noGm_riken$Neg_log_p_val),]


wt_E0771_vs_PY8119_noGm_riken_inf$Neg_log_p_val=303


nrow(wt_E0771_vs_PY8119_noGm_riken)==nrow(wt_E0771_vs_PY8119_noGm_riken_not_inf)+nrow(wt_E0771_vs_PY8119_noGm_riken_inf)

wt_E0771_vs_PY8119_noGm_riken=rbind(wt_E0771_vs_PY8119_noGm_riken_not_inf,wt_E0771_vs_PY8119_noGm_riken_inf)


wt_E0771_vs_PY8119_sig_up=subset(wt_E0771_vs_PY8119_noGm_riken,wt_E0771_vs_PY8119_noGm_riken$log2FoldChange>1 & wt_E0771_vs_PY8119_noGm_riken$padj < 0.01)
wt_E0771_vs_PY8119_sig_dn=subset(wt_E0771_vs_PY8119_noGm_riken,wt_E0771_vs_PY8119_noGm_riken$log2FoldChange < -1 & wt_E0771_vs_PY8119_noGm_riken$padj < 0.01)

wt_E0771_vs_PY8119_sig_up$wt_E0771_vs_PY8119_Direction=ifelse(wt_E0771_vs_PY8119_sig_up$log2FoldChange > 0,"Up")
wt_E0771_vs_PY8119_sig_dn$wt_E0771_vs_PY8119_Direction=ifelse(wt_E0771_vs_PY8119_sig_dn$log2FoldChange < 0,"Down")


wt_E0771_vs_PY8119_sig_up=wt_E0771_vs_PY8119_sig_up[order(- wt_E0771_vs_PY8119_sig_up$log2FoldChange),]
wt_E0771_vs_PY8119_sig_dn=wt_E0771_vs_PY8119_sig_dn[order(wt_E0771_vs_PY8119_sig_dn$log2FoldChange),]

top_15_up=wt_E0771_vs_PY8119_sig_up[1:10,]
top_15_dn=wt_E0771_vs_PY8119_sig_dn[1:10,]

top_15_up$Direction="Up"
top_15_dn$Direction="Down"



top_15_up_dn=rbind(top_15_up,top_15_dn)

#top_15_up_dn=subset(top_15_up_dn,top_15_up_dn$Symbol != "NULL")

### Plotting the top 10 genes as in Ninni's paper

plot_to_save=ggplot(top_15_up_dn, aes(y = reorder(Symbol, log2FoldChange), x = log2FoldChange,size=Neg_log_p_val)) +
  geom_point(data = top_15_dn, aes(y = reorder(Symbol, log2FoldChange), x = log2FoldChange,color=Direction)) +
 geom_point(data = top_15_up, aes(y = reorder(Symbol, log2FoldChange), x = log2FoldChange,color=Direction)) +
 scale_color_manual(values = c("Up" = "red", "Down" = "blue"))+ xlim(-12,25) + geom_vline(xintercept = 0,linetype = "longdash") + ylab("Gene ID") +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle("E0771, WT vs PY8119, WT") + theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=10),axis.title = element_text(size=10))+guides(colour=FALSE)+theme(plot.title = element_text(size=10))+theme(legend.title = element_text(size=10),legend.text = element_text(size=10))


plot_to_save

tiff(file="E0771_vs_PY8119_WT_genes_top_10_genes_by_logFC_plot.tiff",res=300,height = 1500,width = 1500)
grid.draw(plot_to_save)
dev.off()






View(wt_E0771_vs_PY8119_new)


library(ggrepel)
library(grid)


p=ggplot(wt_E0771_vs_PY8119_noGm_riken, aes(baseMean_log, log2FoldChange)) +
  theme_classic(base_size = 16)+
  geom_point(data=wt_E0771_vs_PY8119_noGm_riken, aes(x=baseMean_log, y=log2FoldChange), colour="grey", size=2)
p1 <- p +  geom_point(data = wt_E0771_vs_PY8119_sig_up, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="red")
p2 <- p1 +  geom_point(data = wt_E0771_vs_PY8119_sig_dn, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="blue")
maplot=p2+ggtitle("MA plot for E0771, WT vs PY8119, WT")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title=element_text(size=10))+annotate(geom="text", x=14, y=-15, label="1947 genes up E0771, WT vs PY8119, WT",color="red",size=5)+annotate(geom="text", x=14, y=-17, label="2175 genes down E0771, WT vs PY8119, WT",color="blue",size=5)+xlab("Log2(Mean+1)")+ylim(-20,25)+geom_text_repel(data=top_15_up_dn,aes(x=baseMean_log, y=log2FoldChange,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(0,18)


tiff(file="E0771_vsPY8119_WT_MA_plot.tiff",res=300,height = 3000,width = 3000)
grid.draw(maplot)
dev.off()



p=ggplot(wt_E0771_vs_PY8119_noGm_riken, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 20)+
  geom_point(data=wt_E0771_vs_PY8119_noGm_riken, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = wt_E0771_vs_PY8119_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = wt_E0771_vs_PY8119_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")
volcano=p2+ggtitle("Volcano plot for E0771, WT vs PY8119, WT")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title=element_text(size=14))+annotate(geom="text", x=19.5, y=280, label="1947 genes up",color="red",size=5)+annotate(geom="text", x=19, y=300, label="2175 genes down",color="blue",size=5)+xlab("Log2FC")+ylim(0,320)+geom_text_repel(data=top_15_up_dn,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol,size=5),color="black",arrow=arrow(ends="last",type="open"))+xlim(-20,25)+ylab("-Log10(adj. pval)")

volcano

tiff(file="E0771_vsPY8119_WT_Volcano_plot.tiff",res=300,height = 2000,width = 3000)
grid.draw(volcano)
dev.off()


#boxplot(wt_E0771_vs_PY8119_sig_up_wt_E0771_vs_PY8119_up$log2FoldChange,wt_E0771_vs_PY8119_sig_up_wt_E0771_vs_PY8119_no_change$log2FoldChange,outline=FALSE)

#dn_tab=wt_E0771_vs_PY8119_sig_dn

#top_15_dn=dn_tab[1:15,]

#top_15_up_dn=rbind(top_15_up,top_15_dn)



#rownames(norm_counts)=norm_counts$X
#norm_counts=norm_counts[,-c(1)]
E0771_WT=norm_counts[,c(paste("EOCWM",1:8,sep=""))]
PY8119_WT=norm_counts[,c(paste("PYCWM",1:8,sep=""))]

avg_e0771=apply(E0771_WT,1,mean)
avg_py8119=apply(PY8119_WT,1,mean)

avg_e0771=as.data.frame(avg_e0771)
avg_py8119=as.data.frame(avg_py8119)



rownames(wt_E0771_vs_PY8119_noGm_riken)=wt_E0771_vs_PY8119_noGm_riken$Ensembl

wt_E0771_vs_PY8119_with_exp=cbind(wt_E0771_vs_PY8119_noGm_riken,avg_e0771[rownames(wt_E0771_vs_PY8119_noGm_riken),])
wt_E0771_vs_PY8119_with_exp=cbind(wt_E0771_vs_PY8119_with_exp,avg_py8119[rownames(wt_E0771_vs_PY8119_with_exp),])
colnames(wt_E0771_vs_PY8119_with_exp)[c(14:15)]=c("avg_e0771","avg_py8119")

wt_E0771_vs_PY8119_with_exp$e0771_log=log2(wt_E0771_vs_PY8119_with_exp$avg_e0771+1)
wt_E0771_vs_PY8119_with_exp$py8119_log=log2(wt_E0771_vs_PY8119_with_exp$avg_py8119+1)

rownames(wt_E0771_vs_PY8119_sig_up)=wt_E0771_vs_PY8119_sig_up$Ensembl
rownames(wt_E0771_vs_PY8119_sig_dn)=wt_E0771_vs_PY8119_sig_dn$Ensembl


wt_E0771_vs_PY8119_up_exp=cbind(wt_E0771_vs_PY8119_sig_up,avg_e0771[rownames(wt_E0771_vs_PY8119_sig_up),])
wt_E0771_vs_PY8119_up_exp=cbind(wt_E0771_vs_PY8119_up_exp,avg_py8119[rownames(wt_E0771_vs_PY8119_sig_up),])


wt_E0771_vs_PY8119_dn_exp=cbind(wt_E0771_vs_PY8119_sig_dn,avg_e0771[rownames(wt_E0771_vs_PY8119_sig_dn),])
wt_E0771_vs_PY8119_dn_exp=cbind(wt_E0771_vs_PY8119_dn_exp,avg_py8119[rownames(wt_E0771_vs_PY8119_sig_dn),])

colnames(wt_E0771_vs_PY8119_up_exp)[c(14:15)]=c("avg_e0771","avg_py8119")
colnames(wt_E0771_vs_PY8119_dn_exp)[c(14:15)]=c("avg_e0771","avg_py8119")


wt_E0771_vs_PY8119_up_exp$e0771_log=log2(wt_E0771_vs_PY8119_up_exp$avg_e0771+1)
wt_E0771_vs_PY8119_up_exp$py8119_log=log2(wt_E0771_vs_PY8119_up_exp$avg_py8119+1)

wt_E0771_vs_PY8119_dn_exp$e0771_log=log2(wt_E0771_vs_PY8119_dn_exp$avg_e0771+1)
wt_E0771_vs_PY8119_dn_exp$py8119_log=log2(wt_E0771_vs_PY8119_dn_exp$avg_py8119+1)


rownames(top_15_up_dn)=top_15_up_dn$Ensembl

top_15_up_dn_exp=cbind(top_15_up_dn,avg_e0771[rownames(top_15_up_dn),])
top_15_up_dn_exp=cbind(top_15_up_dn_exp,avg_py8119[rownames(top_15_up_dn_exp),])
colnames(top_15_up_dn_exp)[c(14:15)]=c("avg_e0771","avg_py8119")

top_15_up_dn_exp$avg_e0771_log=log2(top_15_up_dn_exp$avg_e0771+1)
top_15_up_dn_exp$avg_py8119_log=log2(top_15_up_dn_exp$avg_py8119+1)


p=ggplot(wt_E0771_vs_PY8119_with_exp, aes(py8119_log, e0771_log)) +
  theme_classic(base_size = 16)+
  geom_point(data=wt_E0771_vs_PY8119_with_exp, aes(x=py8119_log, y=e0771_log), colour="grey", size=2)
p1 <- p +  geom_point(data = wt_E0771_vs_PY8119_up_exp, aes(x=py8119_log, y=e0771_log) ,size=3,color="red")
p2 <- p1 +  geom_point(data = wt_E0771_vs_PY8119_dn_exp, aes(x=py8119_log, y=e0771_log) ,size=3,color="blue")
correlation_plot=p2+ggtitle("Correlation plot plot for E0771, WT vs PY8119, WT")+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title=element_text(size=10))+annotate(geom="text", x=2, y=13, label="1947 genes up E0771, WT vs PY8119, WT",color="red",size=4)+annotate(geom="text", x=2, y=12, label="2175 genes down E0771, WT vs PY8119, WT",color="blue",size=4)+xlab("Log2(Mean+1) PY8119 WT")+ylim(0,14)+geom_text_repel(data=top_15_up_dn_exp,aes(x=avg_py8119_log,y= avg_e0771_log,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(0,14)+ylim(0,14)+ylab("Log2(Mean+1) E0771 WT")+geom_abline(slope=1,intercept = 0,linetype='dashed')

correlation_plot

tiff(file="E0771_vsPY8119_WT_Correlation_plot.tiff",res=300,height = 3000,width = 3500)
grid.draw(correlation_plot)
dev.off()



