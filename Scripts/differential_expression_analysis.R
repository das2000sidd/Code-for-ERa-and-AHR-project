## This is a prototype of the script I used for the differential expression analysis
## comparing various 6 hour ligand treated MCF-7 cells to 6 hour DMSO treated cells 


DMSO_r1=read.table(file="1-2MCF-DMSO_S7_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
DMSO_r2=read.table(file="2-3MCF-DMSO_S10_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
DMSO_r3=read.table(file="3-4MCF-DMSO_S13_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)


DIM_r1=read.table(file="7-9MCF-DIM_S3_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
DIM_r2=read.table(file="8-11MCF-DIM_S5_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
DIM_r3=read.table(file="9-12MCF-DIM_S8_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)


e2_1=read.table(file="13-17MCF-E2_S18_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
e2_2=read.table(file="14-18MCF-E2_S2_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
e2_3=read.table(file="15-19MCF-E2_S4_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)


res_1=read.table(file="10-13MCF-RES_S11_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
res_2=read.table(file="11-14MCF-RES_S14_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
res_3=read.table(file="12-16MCF-RES_S16_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)


tcdd_1=read.table(file="4-5MCF-TCDD_S15_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
tcdd_2=read.table(file="5-6MCF-TCDD_S17_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
tcdd_3=read.table(file="6-7MCF-TCDD_S1_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)

e2_plus_tcdd_1=read.table(file="16-21MCF-E2-T_S6_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
e2_plus_tcdd_2=read.table(file="17-22MCF-E2-T_S9_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
e2_plus_tcdd_3=read.table(file="18-23MCF-E2-T_S12_R1_001.featureCounts.counts.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)


## It is absolutely critical that the columns of the count matrix and the rows of the column data 
#(information about samples) are in the same order. DESeq2 will not make guesses as to which 
#column of the count matrix belongs to which row of the column data, these must be provided to 
#DESeq2 already in consistent order.
### Due to skewness in darta very difficult to visulaise raw counts
#hist(DMSO_r1$V2, col=rgb(1,0,0,0.5), main="Overlapping Histogram", xlab="Variable",ylim = c(0,1000000))
#hist(e2_plus_tcdd_1$V2, col=rgb(0,0,1,0.5), add=T)

all_counts_combined=cbind(DMSO_r1,DMSO_r2[rownames(DMSO_r1),6],
                          DMSO_r3[rownames(DMSO_r1),6],
                          DIM_r1[rownames(DMSO_r1),6],
                          DIM_r2[rownames(DMSO_r1),6],
                          DIM_r3[rownames(DMSO_r1),6],
                          e2_1[rownames(DMSO_r1),6],
                          e2_2[rownames(DMSO_r1),6],
                          e2_3[rownames(DMSO_r1),6],
                          res_1[rownames(DMSO_r1),6],
                          res_2[rownames(DMSO_r1),6],
                          res_3[rownames(DMSO_r1),6],
                          tcdd_1[rownames(DMSO_r1),6],
                          tcdd_2[rownames(DMSO_r1),6],
                          tcdd_3[rownames(DMSO_r1),6],
                        e2_plus_tcdd_1[rownames(DMSO_r1),6],
                         e2_plus_tcdd_2[rownames(DMSO_r1),6],
                          e2_plus_tcdd_3[rownames(DMSO_r1),6])
#
colnames(all_counts_combined)[6:23]=c("DMSO_r1","DMSO_r2","DMSO_r3",
                                "DIM_r1","DIM_r2","DIM_r3",
                                "e2_1","e2_2","e2_3",
                                "res_1","res_2","res_3",
                                "tcdd_1","tcdd_2","tcdd_3",
                                "e2_plus_tcdd_1","e2_plus_tcdd_2","e2_plus_tcdd_3")
#all_counts_combined=all_counts_combined[1:60605,]

write.table(all_counts_combined[,c(6:23)],file="RNASeq_counts_file_upload_to_GEO.txt",col.names = T,row.names = T,sep="\t",quote = F)

condition_df=as.data.frame(c(rep("DMSO",3),rep("DIM",3),rep("E2",3),rep("RES",3),rep("TCDD",3),rep("E2_plus_TCDD",3)))

colnames(condition_df)="condition"
condition_df$condition=as.factor(condition_df$condition)



library(DESeq2)
library(edgeR)
dds = DESeqDataSetFromMatrix(countData = all_counts_combined[,c(6:23)],
                             colData = condition_df,
                             design = ~ condition)
idx <- rowSums( cpm(dds)) >= 2 ## This is not necessary

head(all_counts_combined[,c(6:8)])



## Exploratory analysis and visualisation

nrow(dds)


vsd <- rlog(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)


library("dplyr")
library("ggplot2")


## Estimating scaling factor per sample
dds <- estimateSizeFactors(dds)


sampleDists <- dist(t(assay(vsd)))
sampleDists


library("pheatmap")
library("RColorBrewer")


sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)




library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)




###PCA plot***
plotPCA(vsd, intgroup = c("condition")) 



###MDS plot***
  mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")  


## MDS plot using the VST data
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color =condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")



library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
top_20_genes=assay(ntd)[select,]
colnames(df)="Condition"
rownames(df)=colnames(top_20_genes)
pheatmap(top_20_genes, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


library(Glimma)
## Running DESeq to fit NB linear models for each gene
dds <- DESeq(dds)

## Visualising PCA plot for batch effect or other discrepancies and seeing how well replicates cluster together
plotMDS(dds)

## Building the results table


res <- results(dds)
res

## DEG for DIM vs DMSO
res_DIM <- results(dds, contrast=c("condition","DIM","DMSO"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_DIM, use.names = TRUE)
res_DIM_df=as.data.frame(res_DIM)
res_DIM_df$Ensembl=rownames(res_DIM_df)
summary(res_DIM)


## DEG for E2+TCDD vs DMSO
res_E2_plus_TCDD <- results(dds, contrast=c("condition","E2_plus_TCDD","DMSO"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_E2_plus_TCDD, use.names = TRUE)
res_E2_plus_TCDD_df=as.data.frame(res_E2_plus_TCDD)
res_E2_plus_TCDD_df$Ensembl=rownames(res_E2_plus_TCDD_df)
summary(res_E2_plus_TCDD)

## DEG for E2 vs DMSO
res_E2 <- results(dds, contrast=c("condition","E2","DMSO"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_E2, use.names = TRUE)
res_E2_df=as.data.frame(res_E2)
res_E2_df$Ensembl=rownames(res_E2_df)
summary(res_E2)

## DEG for RES vs DMSO
res_RES <- results(dds, contrast=c("condition","RES","DMSO"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_RES, use.names = TRUE)
res_RES_df=as.data.frame(res_RES)
res_RES_df$Ensembl=rownames(res_RES_df)
summary(res_RES)

## DEG for TCDD vs DMSO
res_TCDD <- results(dds, contrast=c("condition","TCDD","DMSO"),pAdjustMethod = "BH",format = "DataFrame")
mcols(res_TCDD, use.names = TRUE)
res_TCDD_df=as.data.frame(res_TCDD)
res_TCDD_df$Ensembl=rownames(res_TCDD_df)
summary(res_TCDD)




library(biomaRt)


normalised_counts=counts(dds, normalized=T)

## Pulling out upregulated DEG's using adj.p cutoff of 0.01 and logFC > 0
RES_significant_up=subset(res_RES_df,res_RES_df$padj < 0.01 & res_RES_df$log2FoldChange > 1)
DIM_significant_up=subset(res_DIM_df,res_DIM_df$padj < 0.01 & res_DIM_df$log2FoldChange > 1)
E2_significant_up=subset(res_E2_df,res_E2_df$padj < 0.01 & res_E2_df$log2FoldChange > 1)
TCDD_significant_up=subset(res_TCDD_df,res_TCDD_df$padj < 0.01 & res_TCDD_df$log2FoldChange > 1)
E2_TCDD_significant_up=subset(res_E2_plus_TCDD_df,res_E2_plus_TCDD_df$padj < 0.01 & res_E2_plus_TCDD_df$log2FoldChange > 1)


## Pulling out downregulated DEG's using adj.p cutoff of 0.01 and logFC < -1
RES_significant_dn=subset(res_RES_df,res_RES_df$padj < 0.01 & res_RES_df$log2FoldChange < -1)
DIM_significant_dn=subset(res_DIM_df,res_DIM_df$padj < 0.01 & res_DIM_df$log2FoldChange < -1)
E2_significant_dn=subset(res_E2_df,res_E2_df$padj < 0.01 & res_E2_df$log2FoldChange < -1)
TCDD_significant_dn=subset(res_TCDD_df,res_TCDD_df$padj < 0.01 & res_TCDD_df$log2FoldChange < -1)
E2_TCDD_significant_dn=subset(res_E2_plus_TCDD_df,res_E2_plus_TCDD_df$padj < 0.01 & res_E2_plus_TCDD_df$log2FoldChange < -1)




write.table(res_DIM_df,file="DMSO_vs_DIM_DESEq2_Diff_exp_analysis_Feature_Counts.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(normalised_counts,file="ALL_SAMPLES_NORMALISED_COUNTS.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(res_E2_plus_TCDD_df,file="DMSO_vs_E2_plus_TCDD_DESEq2_Diff_exp_analysis_Feature_Counts.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(res_E2_df,file="DMSO_vs_E2_DESEq2_Diff_exp_analysis_Feature_Counts.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(res_RES_df,file="DMSO_vs_RES_DESEq2_Diff_exp_analysis_Feature_Counts.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(res_TCDD_df,file="DMSO_vs_TCDD_DESEq2_Diff_exp_analysis_Feature_Counts.txt",col.names = T,row.names = T,sep="\t",quote = F)




















