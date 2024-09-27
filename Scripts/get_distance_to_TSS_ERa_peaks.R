
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






