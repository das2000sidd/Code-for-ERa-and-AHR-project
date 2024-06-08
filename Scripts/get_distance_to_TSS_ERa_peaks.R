
### This is a protytype of script I used for plotting the distance to TSS
## for various ChIP seq ligand induced peaks


e2_plus_tcdd=read.table(file="TCDD_plus_E2_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak",header = F,sep="\t",stringsAsFactors = F)
tcdd=read.table(file="TCDD_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak",header = F,sep="\t",stringsAsFactors = F)
dim=read.table(file="DIM_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak",header = F,sep="\t",stringsAsFactors = F)
res=read.table(file="RES_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak",header = F,sep="\t",stringsAsFactors = F)
e2=read.table(file="E2_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak",header = F,sep="\t",stringsAsFactors = F)




colnames(e2_plus_tcdd)[1:3]=c("seqnames","start","end")
colnames(tcdd)[1:3]=c("seqnames","start","end")
colnames(dim)[1:3]=c("seqnames","start","end")
colnames(res)[1:3]=c("seqnames","start","end")
colnames(e2)[1:3]=c("seqnames","start","end")


library("GenomicDistributionsData")
library(GenomicDistributions)
library(GenomicRanges)
library(ChIPseeker)

queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
query = rtracklayer::import(queryFile)

#e2_plus_tcdd$seqnames=paste("chr",e2_plus_tcdd$seqnames,sep="")
#tcdd$seqnames=paste("chr",tcdd$seqnames,sep="")
#dim$seqnames=paste("chr",dim$seqnames,sep="")
#res$seqnames=paste("chr",res$seqnames,sep="")
#e2$seqnames=paste("chr",e2$seqnames,sep="")

e2_plus_tcdd=toGRanges(e2_plus_tcdd[,c(1:3)],format="BED",header=TRUE)
tcdd=toGRanges(tcdd[,c(1:3)],format="BED",header=TRUE)
dim=toGRanges(dim[,c(1:3)],format="BED",header=TRUE)
res=toGRanges(res[,c(1:3)],format="BED",header=TRUE)
e2=toGRanges(e2[,c(1:3)],format="BED",header=TRUE)



e2_TSSdist = calcFeatureDistRefTSS(e2, "hg38")
tcdd_TSSdist = calcFeatureDistRefTSS(tcdd, "hg38")
dim_TSSdist = calcFeatureDistRefTSS(dim, "hg38")
e2_plus_tcdd_TSSdist = calcFeatureDistRefTSS(e2_plus_tcdd, "hg38")
res_TSSdist = calcFeatureDistRefTSS(res, "hg38")



library(ggplot2)

plotFeatureDist(e2_TSSdist,numbers = TRUE)+ggtitle("")+theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1,size=15),axis.text.x = element_text(size=15),axis.title = element_text(size=15))+ylim(0,2000)+ylab("Frequency of peaks")+geom_bar(stat="identity",fill="grey")
plotFeatureDist(tcdd_TSSdist,numbers = TRUE)+ggtitle("")+theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1,size=15),axis.text.x = element_text(size=15),axis.title = element_text(size=15))+ylim(0,2000)+ylab("Frequency of peaks")+geom_bar(stat="identity",fill="grey")
plotFeatureDist(dim_TSSdist,numbers = TRUE)+ggtitle("")+theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1,size=15),axis.text.x = element_text(size=15),axis.title = element_text(size=15))+ylim(0,2000)+ylab("Frequency of peaks")+geom_bar(stat="identity",fill="grey")
plotFeatureDist(e2_plus_tcdd_TSSdist,numbers = TRUE)+ggtitle("")+theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1,size=15),axis.text.x = element_text(size=15),axis.title = element_text(size=15))+ylim(0,2000)+ylab("Frequency of peaks")+geom_bar(stat="identity",fill="grey")
plotFeatureDist(res_TSSdist,numbers = TRUE)+ggtitle("")+theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1,size=15),axis.text.x = element_text(size=15),axis.title = element_text(size=15))+ylim(0,2000)+ylab("Frequency of peaks")+geom_bar(stat="identity",fill="grey")


tiff("Distance to TSS for E2 treated ERa peaks no scale.tiff",width=800)
plotFeatureDist(e2_TSSdist,numbers = TRUE,nbin=200)+ggtitle("")+theme(text = element_text(size=6),axis.text.y = element_blank(),axis.text.x = element_text(size=10),axis.title = element_text(size=6))+ylim(0,400)+ylab("")+geom_bar(stat="identity",fill="grey")
dev.off()

tiff("Distance to TSS for DIM treated ERa peaks no scale.tiff",width=800)
plotFeatureDist(dim_TSSdist,numbers = TRUE,nbin=200)+ggtitle("")+theme(text = element_text(size=6),axis.text.y = element_blank(),axis.text.x = element_text(size=6),axis.title = element_text(size=6))+ylim(0,400)+ylab("")+geom_bar(stat="identity",fill="grey")
dev.off()

tiff("Distance to TSS for RES treated ERa peaks no scale.tiff",width=800)
plotFeatureDist(res_TSSdist,numbers = TRUE,nbin=200)+ggtitle("")+theme(text = element_text(size=6),axis.text.y = element_blank(),axis.text.x = element_text(size=6),axis.title = element_text(size=6))+ylim(0,400)+ylab("")+geom_bar(stat="identity",fill="grey")
dev.off()


tiff("Distance to TSS for E2 treated ERa peaks smaller font.tiff",width=800,height=800,res=500)
plotFeatureDist(e2_TSSdist,numbers = TRUE,nbin=200)+ggtitle("")+theme(text = element_text(size=6),axis.text.y = element_text(angle=90, hjust=1,size=6),axis.text.x = element_text(size=6),axis.title = element_text(size=6))+ylim(0,400)+ylab("Frequency of peaks")+geom_bar(stat="identity",fill="grey")
dev.off()






tiff("Distance to TSS for RES treated ERa peaks.tiff",width=800)
plotFeatureDist(res_TSSdist,numbers = TRUE,nbin=300)+ggtitle("")+theme(text = element_text(size=6),axis.text.y = element_text(angle=90, hjust=1,size=6),axis.text.x = element_text(size=6),axis.title = element_text(size=6))+ylim(0,400)+ylab("Frequency of peaks")+geom_bar(stat="identity",fill="grey")
dev.off()




### overlap each of the regions and find their distances from TSS plot

dim_tcdd_overlap=findOverlapsOfPeaks(dim_grange,tcdd_grange)
dim_tcdd_overlap$venn_cnt


library(ChIPpeakAnno)

makeVennDiagram(list(dim_grange, tcdd_grange), NameOfPeaks=c("DIM peaks, AHR", "TCDD peaks, AHR"),scaled=FALSE, euler.d=FALSE, 
                fill=c("#009E73", "#00FF00"), # circle fill color
                col=c("#D55E00", "#00FF00"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),category.names=c("DIM peaks, AHR", "TCDD peaks, AHR"),cat.cex=1.5,cex=1.5,)


overlap=dim_tcdd_overlap$mergedPeaks
unique=dim_tcdd_overlap$uniquePeaks
unique=as.data.frame(unique)
unique_dim=unique[grep("dim",rownames(unique)),]
unique_tcdd=unique[grep("tcdd",rownames(unique)),]

unique_dim=toGRanges(unique_dim)
unique_tcdd=toGRanges(unique_tcdd)

## unique_dim,unique_tcdd,overlap
TSSdist.overlap = calcFeatureDistRefTSS(overlap, "hg38")
TSSdist.tcdd = calcFeatureDistRefTSS(unique_tcdd, "hg38")
TSSdist.dim = calcFeatureDistRefTSS(unique_dim, "hg38")


library(ggplot2)

plotFeatureDist(TSSdist.overlap, featureName="TSS",numbers = FALSE,size = 1000000)+ggtitle("Distance to TSS for AHR regions overlapping between TCDD and DIM")+ylim(0,20)
plotFeatureDist(TSSdist.tcdd, featureName="TSS",numbers = FALSE,size = 1000000)+ggtitle("Distance to TSS for AHR regions TCDD only")+ylim(0,20)
plotFeatureDist(TSSdist.dim, featureName="TSS",numbers = FALSE,size = 1000000)+ggtitle("Distance to TSS for AHR regions DIM only")+ylim(0,20)




library(AnnotationDbi)
txdb=loadDb("hg38_transcription_txdb.Rdata")
x=as.data.frame(dim_grange)
y=as.data.frame(e2_grange)
z=as.data.frame(res_grange)




dim_tcdd_overlap=findOverlapsOfPeaks(dim_grange,tcdd_grange)
overlap=dim_tcdd_overlap$mergedPeaks
unique_1=dim_tcdd_overlap$uniquePeaks
unique_1=as.data.frame(unique_1)
unique_dim=unique[grep("dim",rownames(unique_1)),]
unique_tcdd=unique[grep("tcdd",rownames(unique_1)),]

unique_dim=makeGRangesFromDataFrame(unique_dim,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames"),
                                    start.field="start",
                                    end.field=c("end", "stop"),
                                    starts.in.df.are.0based=FALSE)


unique_tcdd=makeGRangesFromDataFrame(unique_tcdd,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames"),
                                    start.field="start",
                                    end.field=c("end", "stop"),
                                    starts.in.df.are.0based=FALSE)

overlap_dim_tcdd=makeGRangesFromDataFrame(overlap,
                                    keep.extra.columns=FALSE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames"),
                                    start.field="start",
                                    end.field=c("end", "stop"),
                                    starts.in.df.are.0based=FALSE)






peakAnno.overlap <- annotatePeak(overlap, tssRegion=c(-2000, 2000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno.dim <- annotatePeak(unique_dim, tssRegion=c(-2000, 2000),
                            TxDb=txdb, annoDb="org.Hs.eg.db")


peakAnno.tcdd <- annotatePeak(unique_tcdd, tssRegion=c(-2000, 2000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")


peakAnno.ERa.dim <- annotatePeak(dim_ERa_grange, tssRegion=c(-2000, 2000),
                              TxDb=txdb, annoDb="org.Hs.eg.db")



peakAnno.overlap=as.data.frame(peakAnno.overlap)
peakAnno.dim=as.data.frame(peakAnno.dim)
peakAnno.tcdd=as.data.frame(peakAnno.tcdd)
peakAnno.ERa.dim=as.data.frame(peakAnno.ERa.dim)


peakAnno.overlap.prom=peakAnno.overlap[grep("Promoter",peakAnno.overlap$annotation),]
peakAnno.dim.prom=peakAnno.dim[grep("Promoter",peakAnno.dim$annotation),]
peakAnno.tcdd.prom=peakAnno.tcdd[grep("Promoter",peakAnno.tcdd$annotation),]
peakAnno.ERa.dim.prom=peakAnno.ERa.dim[grep("Promoter",peakAnno.ERa.dim$annotation),]


nrow(peakAnno.tcdd.prom)*100/nrow(peakAnno.tcdd)
nrow(peakAnno.dim.prom)*100/nrow(peakAnno.dim)
nrow(peakAnno.overlap.prom)*100/nrow(peakAnno.overlap)
nrow(peakAnno.ERa.dim.prom)*100/nrow(peakAnno.ERa.dim)



length(unique(peakAnno.tcdd$geneId))
length(unique(peakAnno.dim$geneId))
length(unique(peakAnno.overlap$geneId))
length(unique(peakAnno.ERa.dim$geneId))


rna_dim=read.table(file="DMSO_vs_DIM_DESEq2_Diff_exp_analysis.txt",header = T,sep="\t",stringsAsFactors = F)
rna_e2=read.csv(file="DMSO_vs_E2_DESEq2_Diff_exp_analysis.csv",header = T,stringsAsFactors = F)
rna_tcdd=read.table(file="DMSO_vs_TCDD_DESEq2_Diff_exp_analysis.txt",header = T,sep="\t",stringsAsFactors = F)


rna_dim_sig=subset(rna_dim,rna_dim$padj < 0.01)
rna_e2_sig=subset(rna_e2,rna_e2$padj < 0.01)
rna_tcdd_sig=subset(rna_tcdd,rna_tcdd$padj < 0.01)


length(intersect(peakAnno.dim$geneId,rna_dim_sig$Ensembl)) ## 465
length(intersect(peakAnno.tcdd$geneId,rna_dim_sig$Ensembl)) ## 82
length(intersect(peakAnno.dim$geneId,rna_e2_sig$Ensembl)) ## 652
length(intersect(peakAnno.ERa.dim$geneId,rna_e2_sig$Ensembl)) ## 821
length(intersect(peakAnno.ERa.dim$geneId,rna_dim_sig$Ensembl)) ## 563

Ahr.closest.dim.rnaseq.sig=intersect(peakAnno.dim$geneId,rna_dim_sig$Ensembl)
Er.closest.dim.rnaseq.sig=intersect(peakAnno.ERa.dim$geneId,rna_dim_sig$Ensembl)

length(intersect(Ahr.closest.dim.rnaseq.sig,Er.closest.dim.rnaseq.sig)) ## 334 co regulated by AHR and ER
co.reg = intersect(Ahr.closest.dim.rnaseq.sig,Er.closest.dim.rnaseq.sig)
Ahr.closest.dim.rnaseq.sig.ahr.only= Ahr.closest.dim.rnaseq.sig[Ahr.closest.dim.rnaseq.sig %ni% co.reg]
Er.closest.dim.rnaseq.sig.er.only= Er.closest.dim.rnaseq.sig[Er.closest.dim.rnaseq.sig %ni% co.reg]

length(intersect(Ahr.closest.dim.rnaseq.sig.ahr.only,co.reg))
length(intersect(Er.closest.dim.rnaseq.sig.er.only,co.reg))
length(intersect(Er.closest.dim.rnaseq.sig.er.only,Ahr.closest.dim.rnaseq.sig.ahr.only))

rna_dim.Ahr.closest=rna_dim[rna_dim$Ensembl %in% Ahr.closest.dim.rnaseq.sig.ahr.only,]
rna_dim.Er.closest=rna_dim[rna_dim$Ensembl %in% Er.closest.dim.rnaseq.sig.er.only,]
rna_dim.Ahr.and.Er.closest=rna_dim[rna_dim$Ensembl %in% co.reg,]

write.table(rna_dim.Ahr.closest,file="DIM_altered_genes_closest_to_AHR_binding_site_after_DIM.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(rna_dim.Er.closest,file="DIM_altered_genes_closest_to_ER_binding_site_after_DIM.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(rna_dim.Ahr.and.Er.closest,file="DIM_altered_genes_closest_to_AHR_and_ER_binding_site_after_DIM.txt",col.names = T,row.names = F,sep="\t",quote = F)




x=intersect(peakAnno.dim$geneId,rna_dim_sig$Ensembl)
y=intersect(peakAnno.dim$geneId,rna_e2_sig$Ensembl)

`%ni%` = Negate(`%in%`)
x_unique=x[x %ni% y]

length(intersect(x,y)) ## 412 genes changed after dim are also changed after E2 and are closest to AHR binding site



