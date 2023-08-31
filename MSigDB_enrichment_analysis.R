## This is a prototype of a script I used for GO enrichment analysis of
## all DEgs or DEGs closest to a peak binding site for particular ligand

library(dplyr)

tcdd_r=read.table(file="TCDD_vs_DMSO_RNASeq_diff_exp_run_after_batch_correction.txt",header=T,sep="\t",stringsAsFactors = F)
e2_r=read.table(file="E2_vs_DMSO_RNASeq_diff_exp_run_after_batch_correction.txt",header=T,sep="\t",stringsAsFactors = F)
dim_r=read.table(file="DIM_vs_DMSO_RNASeq_diff_exp_run_after_batch_correction.txt",header=T,sep="\t",stringsAsFactors = F)
res_r=read.table(file="RES_vs_DMSO_RNASeq_diff_exp_run_after_batch_correction.txt",header=T,sep="\t",stringsAsFactors = F)
e2_tcdd_r=read.table(file="E2_plus_TCDD_vs_DMSO_RNASeq_diff_exp_run_after_batch_correction.txt",header=T,sep="\t",stringsAsFactors = F)

bakcground_genes_detected=read.table(file="Genes_used_for_diff_exp_run_after_batch_correction.txt",header = F,sep="\t",stringsAsFactors = F)


library(org.Hs.eg.db)

e2_r$Entrez <- mapIds(org.Hs.eg.db, e2_r$Ensembl,keytype="ENSEMBL", column="ENTREZID")
dim_r$Entrez <- mapIds(org.Hs.eg.db, dim_r$Ensembl,keytype="ENSEMBL", column="ENTREZID")
res_r$Entrez <- mapIds(org.Hs.eg.db, res_r$Ensembl,keytype="ENSEMBL", column="ENTREZID")
e2_tcdd_r$Entrez <- mapIds(org.Hs.eg.db, e2_tcdd_r$Ensembl,keytype="ENSEMBL", column="ENTREZID")
tcdd_r$Entrez <- mapIds(org.Hs.eg.db, tcdd_r$Ensembl,keytype="ENSEMBL", column="ENTREZID")
bakcground_genes_detected$Entrez <- mapIds(org.Hs.eg.db, bakcground_genes_detected$V1,keytype="ENSEMBL", column="ENTREZID")


e2_r$Symbol <- mapIds(org.Hs.eg.db, e2_r$Entrez,keytype="ENTREZID", column="SYMBOL")
dim_r$Symbol <- mapIds(org.Hs.eg.db, dim_r$Entrez,keytype="ENTREZID", column="SYMBOL")
res_r$Symbol <- mapIds(org.Hs.eg.db, res_r$Entrez,keytype="ENTREZID", column="SYMBOL")
e2_tcdd_r$Symbol <- mapIds(org.Hs.eg.db, e2_tcdd_r$Entrez,keytype="ENTREZID", column="SYMBOL")
tcdd_r$Symbol <- mapIds(org.Hs.eg.db, tcdd_r$Entrez,keytype="ENTREZID", column="SYMBOL")
bakcground_genes_detected$Symbol <- mapIds(org.Hs.eg.db, bakcground_genes_detected$Entrez,keytype="ENTREZID", column="SYMBOL")



e2_r$Symbol=as.character(e2_r$Symbol)
dim_r$Symbol=as.character(dim_r$Symbol)
res_r$Symbol=as.character(res_r$Symbol)
e2_tcdd_r$Symbol=as.character(e2_tcdd_r$Symbol)
bakcground_genes_detected$Symbol=as.character(bakcground_genes_detected$Symbol)
tcdd_r$Symbol=as.character(tcdd_r$Symbol)


e2_pos=subset(e2_r,e2_r$padj < 0.01 & e2_r$log2FoldChange > 1)
e2_neg=subset(e2_r,e2_r$padj < 0.01 & e2_r$log2FoldChange < -1)
e2_sig=subset(e2_r,e2_r$padj < 0.01 & abs(e2_r$log2FoldChange) > 1)



dim_pos=subset(dim_r,dim_r$padj < 0.01 & dim_r$log2FoldChange > 1)
dim_neg=subset(dim_r,dim_r$padj < 0.01 & dim_r$log2FoldChange < -1)
dim_sig=subset(dim_r,dim_r$padj < 0.01 & abs(dim_r$log2FoldChange) > 1)


res_pos=subset(res_r,res_r$padj < 0.01 & res_r$log2FoldChange > 1)
res_neg=subset(res_r,res_r$padj < 0.01 & res_r$log2FoldChange < -1)
res_sig=subset(res_r,res_r$padj < 0.01 & abs(res_r$log2FoldChange) > 1)


e2_tcdd_pos=subset(e2_tcdd_r,e2_tcdd_r$padj < 0.01 & e2_tcdd_r$log2FoldChange > 1)
e2_tcdd_neg=subset(e2_tcdd_r,e2_tcdd_r$padj < 0.01 & e2_tcdd_r$log2FoldChange < -1)
e2_tcdd_sig=subset(e2_tcdd_r,e2_tcdd_r$padj < 0.01 & abs(e2_tcdd_r$log2FoldChange) > 1)


tcdd_sig=subset(tcdd_r,tcdd_r$padj < 0.01 & abs(tcdd_r$log2FoldChange) > 1)


intersect(e2_sig$Ensembl,e2_tcdd_sig$Ensembl)
intersect(e2_pos$Ensembl,e2_tcdd_pos$Ensembl)
intersect(e2_neg$Ensembl,e2_tcdd_neg$Ensembl)
intersect(dim_sig$Ensembl,e2_tcdd_sig$Ensembl)


e2_pos_genes=unique(e2_pos$Entrez)
dim_pos_genes=unique(dim_pos$Entrez)
res_pos_genes=unique(res_pos$Entrez)

e2_neg_genes=unique(e2_neg$Entrez)
dim_neg_genes=unique(dim_neg$Entrez)
res_neg_genes=unique(res_neg$Entrez)

e2_neg_pos=unique(c(e2_pos_genes,e2_neg_genes))
dim_neg_pos=unique(c(dim_pos_genes,dim_neg_genes))
res_neg_pos=unique(c(res_pos_genes,res_neg_genes))
e2_tcdd_neg_pos=unique(c(e2_tcdd_sig$Entrez))




dim_er=read.table(file="DIM_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
e2_er=read.table(file="E2_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
res_er=read.table(file="RES_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
dim_ahr=read.table(file="DIM_vs_DMSO_AHR_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
e2_TCDD_ahr=read.table(file="TCDD_plus_E2_vs_DMSO_AHR_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
e2_TCDD_er=read.table(file="TCDD_plus_E2_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
tcdd_ahr=read.table(file="TCDD_vs_DMSO_AHR_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)




colnames(dim_er)[1:3]=c("seqnames","start","end")
colnames(e2_er)[1:3]=c("seqnames","start","end")
colnames(res_er)[1:3]=c("seqnames","start","end")
colnames(dim_ahr)[1:3]=c("seqnames","start","end")
colnames(e2_TCDD_ahr)[1:3]=c("seqnames","start","end")
colnames(e2_TCDD_er)[1:3]=c("seqnames","start","end")
colnames(tcdd_ahr)[1:3]=c("seqnames","start","end")





chr_keep=1:22
chr_all=c(chr_keep,"X")


dim_er=dim_er[dim_er$seqnames %in% chr_all,]
e2_er=e2_er[e2_er$seqnames %in% chr_all,]
res_er=res_er[res_er$seqnames %in% chr_all,]
dim_ahr=dim_ahr[dim_ahr$seqnames %in% chr_all,]
e2_TCDD_ahr=e2_TCDD_ahr[e2_TCDD_ahr$seqnames %in% chr_all,]
e2_TCDD_er=e2_TCDD_er[e2_TCDD_er$seqnames %in% chr_all,]
tcdd_ahr=tcdd_ahr[tcdd_ahr$seqnames %in% chr_all,]




library(ChIPseeker)
library(ChIPpeakAnno)
library(AnnotationDbi)

dim_er=toGRanges(dim_er[,c(1:3)],format="BED",header=TRUE)
e2_er=toGRanges(e2_er[,c(1:3)],format="BED",header=TRUE)
res_er=toGRanges(res_er[,c(1:3)],format="BED",header=TRUE)
dim_ahr=toGRanges(dim_ahr[,c(1:3)],format="BED",header=TRUE)
e2_TCDD_ahr=toGRanges(e2_TCDD_ahr[,c(1:3)],format="BED",header=TRUE)
e2_TCDD_er=toGRanges(e2_TCDD_er[,c(1:3)],format="BED",header=TRUE)
tcdd_ahr=toGRanges(tcdd_ahr[,c(1:3)],format="BED",header=TRUE)



txdb=loadDb("hg38_transcription_txdb.Rdata")

peakAnno.dim.er <- annotatePeak(dim_er, tssRegion=c(-1000, 1000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.e2.er <- annotatePeak(e2_er, tssRegion=c(-1000, 1000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.res.er <- annotatePeak(res_er, tssRegion=c(-1000, 1000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.dim.ahr <- annotatePeak(dim_ahr, tssRegion=c(-1000, 1000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.e2.plus.tcdd.ahr <- annotatePeak(e2_TCDD_ahr, tssRegion=c(-1000, 1000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.e2.plus.tcdd.er <- annotatePeak(e2_TCDD_er, tssRegion=c(-1000, 1000),
                                          TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.tcdd.ahr <- annotatePeak(tcdd_ahr, tssRegion=c(-1000, 1000),
                                         TxDb=txdb, annoDb="org.Hs.eg.db")




peakAnno.dim.er=as.data.frame(peakAnno.dim.er)
peakAnno.e2.er=as.data.frame(peakAnno.e2.er)
peakAnno.res.er=as.data.frame(peakAnno.res.er)
peakAnno.dim.ahr=as.data.frame(peakAnno.dim.ahr)
peakAnno.e2.plus.tcdd.ahr=as.data.frame(peakAnno.e2.plus.tcdd.ahr)
peakAnno.e2.plus.tcdd.er=as.data.frame(peakAnno.e2.plus.tcdd.er)
peakAnno.tcdd.ahr=as.data.frame(peakAnno.tcdd.ahr)



dim_up_dn_er=intersect(dim_sig$Ensembl,peakAnno.dim.er$geneId)
e2_up_dn_er=intersect(e2_sig$Ensembl,peakAnno.e2.er$geneId)
res_up_dn_er=intersect(res_sig$Ensembl,peakAnno.res.er$geneId)
dim_up_dn_ahr=intersect(dim_sig$Ensembl,peakAnno.dim.ahr$geneId)
e2_tcdd_up_dn_ahr=intersect(e2_tcdd_sig$Ensembl,peakAnno.e2.plus.tcdd.ahr$geneId)
e2_tcdd_up_dn_er=intersect(e2_tcdd_sig$Ensembl,peakAnno.e2.plus.tcdd.er$geneId)
tcdd_up_dn_ahr=intersect(tcdd_sig$Ensembl,peakAnno.tcdd.ahr$geneId)


e2_up_er=intersect(e2_pos$Ensembl,peakAnno.e2.er$geneId)
res_up_er=intersect(res_pos$Ensembl,peakAnno.res.er$geneId)


e2_dn_er=intersect(e2_neg$Ensembl,peakAnno.e2.er$geneId)
res_dn_er=intersect(res_neg$Ensembl,peakAnno.res.er$geneId)


library(eulerr)

dim_up_dn_er=subset(dim_r,dim_r$Ensembl %in% dim_up_dn_er)
e2_up_dn_er=subset(e2_r,e2_r$Ensembl %in% e2_up_dn_er)
res_up_dn_er=subset(res_r,res_r$Ensembl %in% res_up_dn_er)
e2_tcdd_up_dn_er=subset(e2_tcdd_r,e2_tcdd_r$Ensembl %in% e2_tcdd_up_dn_er)



dim_up_dn_er_ahr = intersect(dim_up_dn_ahr$Ensembl,dim_up_dn_er$Ensembl)

dim_up_dn_er_ahr=subset(dim_r,dim_r$Ensembl %in% dim_up_dn_er_ahr)


e2_up=subset(e2_r,e2_r$padj < 0.01 & e2_r$log2FoldChange > 1)
e2_dn=subset(e2_r,e2_r$padj < 0.01 & e2_r$log2FoldChange < -1)

dim_dn=subset(dim_r,dim_r$padj < 0.01 & dim_r$log2FoldChange < -1)
dim_up=subset(dim_r,dim_r$padj < 0.01 & dim_r$log2FoldChange > 1)


dim_up_dn_er = intersect(dim_sig$Ensembl,peakAnno.dim.er$geneId)
dim_up_dn_ahr = intersect(dim_sig$Ensembl,peakAnno.dim.ahr$geneId)


dim_up_dn_er = intersect(dim_sig$Ensembl,peakAnno.dim.er$geneId)

dim_up_dn_er_ahr=intersect(dim_up_dn_ahr,dim_up_dn_er)


dim_sig_er=subset(dim_sig,dim_sig$Ensembl %in% peakAnno.dim.er$geneId)

dim_up_dn_er_ahr=subset(dim_r,dim_r$Ensembl %in% dim_up_dn_er_ahr)


dim_sig_up_er=subset(dim_pos,dim_pos$Ensembl %in% peakAnno.dim.er$geneId)
dim_sig_dn_er=subset(dim_neg,dim_neg$Ensembl %in% peakAnno.dim.er$geneId)


dim_sig_up_ahr=subset(dim_pos,dim_pos$Ensembl %in% peakAnno.dim.ahr$geneId)
dim_sig_dn_ahr=subset(dim_neg,dim_neg$Ensembl %in% peakAnno.dim.ahr$geneId)



e2_sig_up_er=subset(e2_pos,e2_pos$Ensembl %in% peakAnno.e2.er$geneId)
e2_sig_dn_er=subset(e2_neg,e2_neg$Ensembl %in% peakAnno.e2.er$geneId)
e2_sig_er=subset(e2_sig,e2_sig$Ensembl %in% peakAnno.e2.er$geneId)


res_sig_er=subset(res_sig,res_sig$Ensembl %in% peakAnno.res.er$geneId)


e2_tcdd_sig_er=subset(e2_tcdd_sig,e2_tcdd_sig$Ensembl %in% peakAnno.e2.plus.tcdd.er$geneId)


e2_tcdd_sig_ahr=subset(e2_tcdd_sig,e2_tcdd_sig$Ensembl %in% peakAnno.e2.plus.tcdd.ahr$geneId)

e2_tcdd_sig_er_ahr=intersect(e2_tcdd_sig_ahr$Ensembl,e2_tcdd_sig_er$Ensembl)


dim_sig_ahr=subset(dim_sig,dim_sig$Ensembl %in% peakAnno.dim.ahr$geneId)


dim_sig_er_ahr=intersect(dim_sig_ahr$Ensembl,dim_sig_er$Ensembl)


### enrichment analysis using MSigDB
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(magrittr)

msigdbr_species()

hs_msigdb_df <- msigdbr(species = "Homo sapiens")

head(hs_msigdb_df)

#Filter the human data frame to the KEGG pathways that are included in the
# curated gene sets
hs_GO_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C5", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("GO:BP","GO:CC","GO:MF") # This is because we only want KEGG pathways
  )

## DO NOT USE
hs_KEGG_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("CP:KEGG") # This is because we only want KEGG pathways
  )

## DO NOT USE
hs_hallmark <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "H" # This is to filter only to the C2 curated gene sets
     # This is because we only want KEGG pathways
  )



## ONLY RESULTS TO USE
GO_ora_results_e2 <- enricher(
  gene = e2_sig$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_e2@result)

enrich_plot <- enrichplot::dotplot(GO_ora_results_e2, showCategory=15,font.size=10,title="E2 vs DMSO DEG enrichment using MsigDB top 20 most significant GO terms by adj. p value",orderBy= "p.adjust")
tiff("E2_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()



e2_sig_er=subset(e2_sig,e2_sig$Ensembl %in% peakAnno.e2.er$geneId)


hallmark_ora_results_e2_er <- enricher(
  gene = e2_sig_er$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)


enrich_plot <- enrichplot::dotplot(hallmark_ora_results_e2_er, showCategory=15,font.size=10,title="E2 vs DMSO significant genes closest to ERa enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot
tiff("E2_ERa_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()




dim_sig_er=subset(dim_sig,dim_sig$Ensembl %in% peakAnno.dim.er$geneId)


GO_ora_results_dim <- enricher(
  gene = dim_sig$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_dim@result)

enrich_plot <- enrichplot::dotplot(GO_ora_results_dim, showCategory=15,font.size=10,title="DIM vs DMSO significant genes enrichment using MsigDB top 20 GO terms",orderBy="p.adjust")
enrich_plot
tiff("DIM_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()



GO_ora_results_dim_er <- enricher(
  gene = dim_sig_er$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_dim_er@result)

enrich_plot <- enrichplot::dotplot(GO_ora_results_dim_er, showCategory=15,font.size=10,title="DIM vs DMSO significant genes closest to ERa enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot
tiff("DIM_ERa_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()






GO_ora_results_res <- enricher(
  gene = res_sig$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_res@result)

enrich_plot <- enrichplot::dotplot(GO_ora_results_res, showCategory=15,font.size=10,title="RES vs DMSO significant genes enrichment using MsigDB top 20 GO terms",orderBy="p.adjust")
enrich_plot
tiff("RES_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()


res_sig_er=subset(res_sig,res_sig$Ensembl %in% peakAnno.res.er$geneId)




GO_ora_results_res_er <- enricher(
  gene = res_sig_er$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_res_er@result)

enrich_plot <- enrichplot::dotplot(GO_ora_results_res_er, showCategory=15,font.size=10,title="RES vs DMSO significant genes closest to ERa enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot
tiff("RES_ERa_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()





dim_sig_ahr=subset(dim_sig,dim_sig$Ensembl %in% peakAnno.dim.ahr$geneId)


GO_ora_results_dim_ahr <- enricher(
  gene = dim_sig_ahr$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_dim_ahr@result)

enrich_plot <- enrichplot::dotplot(GO_ora_results_dim_ahr, showCategory=15,font.size=10,title="DIM vs DMSO significant genes closest to AHR enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot
tiff("DIM_AHR_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()







GO_ora_results_e2_tcdd <- enricher(
  gene = e2_tcdd_sig$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_e2_tcdd@result)

enrich_plot <- enrichplot::dotplot(GO_ora_results_e2_tcdd, showCategory=15,font.size=10,title="E2+TCDD vs DMSO significant genes enrichment using MsigDB top 20 GO terms",orderBy="p.adjust")
enrich_plot
tiff("E2_plus_TCDD_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()



e2_plus_tcdd_sig_ahr=subset(e2_tcdd_sig,e2_tcdd_sig$Ensembl %in% peakAnno.e2.plus.tcdd.ahr$geneId)


GO_ora_results_e2_tcdd_ahr <- enricher(
  gene = e2_plus_tcdd_sig_ahr$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

enrich_plot <- enrichplot::dotplot(GO_ora_results_e2_tcdd_ahr, showCategory=15,font.size=10,title="E2+TCDD vs DMSO significant genes closest to AHR enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot
tiff("E2_plus_TCDD_AHR_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()



e2_plus_tcdd_sig_er=subset(e2_tcdd_sig,e2_tcdd_sig$Ensembl %in% peakAnno.e2.plus.tcdd.er$geneId)

e2_plus_tcdd_sig_ahr_er=intersect(e2_plus_tcdd_sig_ahr$Ensembl,e2_plus_tcdd_sig_er$Ensembl)





GO_ora_results_e2_tcdd_ahr_er <- enricher(
  gene = e2_plus_tcdd_sig_ahr_er, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

enrich_plot <- enrichplot::dotplot(GO_ora_results_e2_tcdd_ahr_er, showCategory=15,font.size=10,title="E2+TCDD vs DMSO significant genes closest to AHR and ERa enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot
tiff("E2_plus_TCDD_ERa_AHR_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()





GO_ora_results_e2_tcdd_er <- enricher(
  gene = e2_plus_tcdd_sig_er$Ensembl, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

enrich_plot <- enrichplot::dotplot(GO_ora_results_e2_tcdd_er, showCategory=15,font.size=10,title="E2+TCDD vs DMSO significant genes closest to ERa enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot
tiff("E2_plus_TCDD_ERa_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()





dim_sig_ahr_er=intersect(dim_sig_ahr$Ensembl,dim_sig_er$Ensembl)





GO_ora_results_dim_er_ahr <- enricher(
  gene = dim_sig_ahr_er, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = bakcground_genes_detected$V1,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)

enrich_plot <- enrichplot::dotplot(GO_ora_results_dim_er_ahr, showCategory=15,font.size=10,title="DIM vs DMSO significant genes closest to AHR and ERa enrichment using MsigDB GO terms",orderBy="p.adjust")
enrich_plot

tiff("DIM_ERa_AHR_GO_enrichment_derived_from_MSigDB_database.tiff",width=800)
enrich_plot ## Dis
dev.off()






View(GO_ora_results_e2@result)

write.table(GO_ora_results_e2@result,file="E2_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(hallmark_ora_results_e2_er@result,file="E2_ERa_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_e2_tcdd_er@result,file="E2_plus_TCDD_ERa_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_dim@result,file="DIM_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_dim_er@result,file="DIM_ERa_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_res@result,file="RES_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_res_er@result,file="RES_ERa_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_dim_ahr@result,file="DIM_AHR_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_e2_tcdd@result,file="E2_plus_TCDD_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_e2_tcdd_ahr@result,file="E2_plus_TCDD_AHR_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_e2_tcdd_ahr_er@result,file="E2_plus_TCDD_AHR_ERa_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_e2_tcdd_er@result,file="E2_plus_TCDD_ERa_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(GO_ora_results_dim_er_ahr@result,file="DIM_ERa_AHR_GO_enrichment_derived_from_MSigDB_database_symbol.txt",col.names = T,row.names = T,sep="\t",quote = F)






