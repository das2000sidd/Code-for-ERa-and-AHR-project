## This is a prototype of a script I used for plotting
## upset plot for overlapping differentially expressed genes
## for various ligand treated samples relative to control DMSO

e2_r=read.table(file="DMSO_vs_E2_DESEq2_Diff_exp_analysis.txt",header=T,sep="\t",stringsAsFactors = F)
dim_r=read.table(file="DMSO_vs_DIM_DESEq2_Diff_exp_analysis.txt",header=T,sep="\t",stringsAsFactors = F)
res_r=read.table(file="DMSO_vs_RES_DESEq2_Diff_exp_analysis.txt",header=T,sep="\t",stringsAsFactors = F)


e2_sig=subset(e2_r,e2_r$padj < 0.01 & abs(e2_r$log2FoldChange) > 1)
dim_sig=subset(dim_r,dim_r$padj < 0.01 & abs(dim_r$log2FoldChange) > 1)
res_sig=subset(res_r,res_r$padj < 0.01 & abs(res_r$log2FoldChange) > 1)

e2_up=subset(e2_sig,e2_sig$log2FoldChange > 1)
dim_up=subset(dim_sig,dim_sig$log2FoldChange > 1)
res_up=subset(res_sig,res_sig$log2FoldChange > 1)

e2_dn=subset(e2_sig,e2_sig$log2FoldChange <  -1)
dim_dn=subset(dim_sig,dim_sig$log2FoldChange < -1)
res_dn=subset(res_sig,res_sig$log2FoldChange < -1)





e2_up_genes=e2_up$Ensembl
dim_up_genes=dim_up$Ensembl
res_up_genes=res_up$Ensembl

e2_dn_genes=e2_dn$Ensembl
dim_dn_genes=dim_dn$Ensembl
res_dn_genes=res_dn$Ensembl

all_up_genes=c(e2_up_genes,dim_up_genes,
               res_up_genes)
all_up_genes=unique(all_up_genes)

all_dn_genes=c(e2_dn_genes,dim_dn_genes,
               res_dn_genes)
all_dn_genes=unique(all_dn_genes)

e2_up_genes=as.data.frame(e2_up_genes)
dim_up_genes=as.data.frame(dim_up_genes)
res_up_genes=as.data.frame(res_up_genes)

rownames(e2_up_genes)=e2_up_genes$e2_up_genes
rownames(dim_up_genes)=dim_up_genes$dim_up_genes
rownames(res_up_genes)=res_up_genes$res_up_genes

e2_dn_genes=as.data.frame(e2_dn_genes)
dim_dn_genes=as.data.frame(dim_dn_genes)
res_dn_genes=as.data.frame(res_dn_genes)

rownames(e2_dn_genes)=e2_dn_genes$e2_dn_genes
rownames(dim_dn_genes)=dim_dn_genes$dim_dn_genes
rownames(res_dn_genes)=res_dn_genes$res_dn_genes



e2_up_genes$Type=1
dim_up_genes$Type=1
res_up_genes$Type=1

e2_dn_genes$Type=1
dim_dn_genes$Type=1
res_dn_genes$Type=1




all_up_genes=as.data.frame(all_up_genes)
rownames(all_up_genes)=all_up_genes$all_up_genes

all_dn_genes=as.data.frame(all_dn_genes)
rownames(all_dn_genes)=all_dn_genes$all_dn_genes




all_up_genes_table=cbind(all_up_genes,
                         e2_up_genes[rownames(all_up_genes),2],
                         dim_up_genes[rownames(all_up_genes),2],
                         res_up_genes[rownames(all_up_genes),2])

colnames(all_up_genes_table)[c(2:4)]=c("E2","DIM","RES")

all_up_genes_table[is.na(all_up_genes_table)]=0


all_dn_genes_table=cbind(all_dn_genes,
                         e2_dn_genes[rownames(all_dn_genes),2],
                         dim_dn_genes[rownames(all_dn_genes),2],
                         res_dn_genes[rownames(all_dn_genes),2])

colnames(all_dn_genes_table)[c(2:4)]=c("E2","DIM","RES")


all_dn_genes_table[is.na(all_dn_genes_table)]=0

colnames(all_up_genes_table)[1]="Genes"
colnames(all_dn_genes_table)[1]="Genes"


all_sig_Genes=rbind(all_up_genes_table,all_dn_genes_table)


library(UpSetR)

upset(all_up_genes_table[,c(2:4)],order.by = c("freq"),main.bar.color = c("black"),set_size.show=TRUE,set_size.numbers_size = 16,text.scale = c(3,2,2,1,2,2),set_size.scale_max = 750,point.size = 3,mainbar.y.label = "Intersection Size")


upset(all_dn_genes_table[,c(2:4)],order.by = c("freq"),main.bar.color = c("black"),set_size.show=TRUE,set_size.numbers_size = 16,text.scale = c(3,2,2,1,2,2),set_size.scale_max = 500,point.size = 3,mainbar.y.label = "Intersection Size")


upset(all_sig_Genes[,c(2:4)],order.by = c("freq"),main.bar.color = c("black"),set_size.show=TRUE,set_size.numbers_size = 16,text.scale = c(3,2,2,1,2,2),set_size.scale_max = 1250,point.size = 3,mainbar.y.label = "Intersection Size")



