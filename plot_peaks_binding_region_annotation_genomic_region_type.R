## This is a prototype of a script I used for plotting
## piechart of genomewide annotation of ligand induced peak regions
## as well as overlapping ligand induced peaks nearest gene
## for estrogen receptor alpha ligand induced peaks


library(ChIPseeker)
library(ChIPpeakAnno)
library(AnnotationDbi)

dim=read.table(file="DIM_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
e2=read.table(file="E2_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
tcdd=read.table(file="TCDD_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
e2_tcdd=read.table(file="TCDD_plus_E2_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)
res=read.table(file="RES_vs_DMSO_ERa_macs3_0.0001_peaks.narrowPeak.no.chr.word",header = F,sep="\t",stringsAsFactors = F)


colnames(e2_tcdd)[1:3]=c("seqnames","start","end")
colnames(dim)[1:3]=c("seqnames","start","end")
colnames(res)[1:3]=c("seqnames","start","end")
colnames(e2)[1:3]=c("seqnames","start","end")
colnames(tcdd)[1:3]=c("seqnames","start","end")



chr_keep=1:22
chr_all=c(chr_keep,"X","Y")
#chr_all=paste("chr",chr_all,sep="")

dim=dim[dim$seqnames %in% chr_all,]
e2=e2[e2$seqnames %in% chr_all,]
res=res[res$seqnames %in% chr_all,]
e2_tcdd=e2_tcdd[e2_tcdd$seqnames %in% chr_all,]
tcdd=tcdd[tcdd$seqnames %in% chr_all,]



e2=toGRanges(e2[,c(1:3)],format="BED",header=TRUE)
dim=toGRanges(dim[,c(1:3)],format="BED",header=TRUE)
res=toGRanges(res[,c(1:3)],format="BED",header=TRUE)
e2_tcdd=toGRanges(e2_tcdd[,c(1:3)],format="BED",header=TRUE)
tcdd=toGRanges(tcdd[,c(1:3)],format="BED",header=TRUE)


txdb=loadDb("hg38_transcription_txdb.Rdata")


peakAnno.e2 <- annotatePeak(e2, tssRegion=c(-1000, 1000),
                                   TxDb=txdb, annoDb="org.Hs.eg.db") ## 

peakAnno.dim <- annotatePeak(dim, tssRegion=c(-1000, 1000),
                              TxDb=txdb, annoDb="org.Hs.eg.db") ## 

peakAnno.res <- annotatePeak(res, tssRegion=c(-1000, 1000),
                             TxDb=txdb, annoDb="org.Hs.eg.db") ## 

peakAnno.e2.tcdd <- annotatePeak(e2_tcdd, tssRegion=c(-1000, 1000),
                             TxDb=txdb, annoDb="org.Hs.eg.db") ## 

peakAnno.tcdd <- annotatePeak(tcdd, tssRegion=c(-1000, 1000),
                                 TxDb=txdb, annoDb="org.Hs.eg.db") ## 



tiff("E2_ERa_genomic_region_distribution.tiff",width=800)
plotAnnoPie(peakAnno.e2,main="Distribution of ERα binding sites by genomic region after E2",cex=1)
dev.off()

tiff("DIM_ERa_genomic_region_distribution.tiff",width=800)
plotAnnoPie(peakAnno.dim,main="Distribution of ERα binding sites by genomic region after DIM",cex=1)
dev.off()

tiff("RES_ERa_genomic_region_distribution.tiff",width=800)
plotAnnoPie(peakAnno.res,main="Distribution of ERα binding sites by genomic region after RES",cex=1)
dev.off()




plotAnnoPie(peakAnno.dim,main="Distribution of ERα binding sites by genomic region after DIM")
plotAnnoPie(peakAnno.res,main="Distribution of ERα binding sites by genomic region after RES")
plotAnnoPie(peakAnno.e2.tcdd,main="Distribution of ERα binding sites by genomic region after E2+TCDD")
plotAnnoPie(peakAnno.tcdd,main="Distribution of ERα binding sites by genomic region after TCDD")


peakAnno.e2=as.data.frame(peakAnno.e2)
peakAnno.dim=as.data.frame(peakAnno.dim)
peakAnno.res=as.data.frame(peakAnno.res)
peakAnno.e2.tcdd=as.data.frame(peakAnno.e2.tcdd)
peakAnno.tcdd=as.data.frame(peakAnno.tcdd)

write.table(peakAnno.e2,file="E2_ERa_peaks_annotation.txt",sep="\t",col.names = T,quote = F,row.names = F)
write.table(peakAnno.dim,file="DIM_ERa_peaks_annotation.txt",sep="\t",col.names = T,quote = F,row.names = F)
write.table(peakAnno.res,file="RES_ERa_peaks_annotation.txt",sep="\t",col.names = T,quote = F,row.names = F)



length(unique(peakAnno.e2$geneId))
length(unique(peakAnno.dim$geneId))
length(unique(peakAnno.res$geneId))
length(unique(peakAnno.e2.tcdd$geneId))
length(unique(peakAnno.tcdd$geneId))


peakAnno.e2.prom=peakAnno.e2[grep("Promoter",peakAnno.e2$annotation),] ## 
nrow(peakAnno.e2.prom) ## 7022
peakAnno.dim.prom=peakAnno.dim[grep("Promoter",peakAnno.dim$annotation),]
nrow(peakAnno.dim.prom)  ## 5289
peakAnno.res.prom=peakAnno.res[grep("Promoter",peakAnno.res$annotation),]
nrow(peakAnno.res.prom)  # 5264
peakAnno.e2.tcdd.prom=peakAnno.e2.tcdd[grep("Promoter",peakAnno.e2.tcdd$annotation),]
nrow(peakAnno.e2.tcdd.prom) ## 8435
peakAnno.tcdd.prom=peakAnno.tcdd[grep("Promoter",peakAnno.tcdd$annotation),]
nrow(peakAnno.tcdd.prom) ## 132


peakAnno.e2.prom.genes=unique(peakAnno.e2.prom$geneId)
length(peakAnno.e2.prom.genes) ## 5210
peakAnno.dim.prom.genes=unique(peakAnno.dim.prom$geneId)
length(peakAnno.dim.prom.genes)  ## 4102
peakAnno.res.prom.genes=unique(peakAnno.res.prom$geneId)
length(peakAnno.res.prom.genes)  # 4056
peakAnno.e2.plus.tcdd.genes=unique(peakAnno.e2.tcdd.prom$geneId)
length(peakAnno.e2.plus.tcdd.genes) ## 6076
peakAnno.tcdd.genes=unique(peakAnno.tcdd.prom$geneId)
length(peakAnno.tcdd.genes) ## 119


e2.genes = unique(peakAnno.e2$geneId)
dim.genes = unique(peakAnno.dim$geneId)
res.genes = unique(peakAnno.res$geneId)
e2.tcdd.genes = unique(peakAnno.e2.tcdd$geneId)
tcdd.genes = unique(peakAnno.tcdd$geneId)


library(VennDiagram)
library(eulerr)


e2.dim=venn.diagram(list("E2,ERα"=e2.genes,"DIM,ERα"=dim.genes),cex=2,fill=c("white","white"),main.cex=c(1),cat.cex=c(1.5,1.5),cat.pos=c(-30,45),cat.dist=c(.04,.04),euler.d=TRUE,filename="E2_DIM_ERa_nearest_gene_overlap_high_res.png",imagetype = "png",height=1000,width=1000,resolution = 300)
grid.draw(e2.dim)
e2.dim=venn.diagram(list("E2,ERα"=e2.genes,"DIM,ERα"=dim.genes),cex=2,fill=c("white","white"),main.cex=c(1),cat.cex=c(1.5,1.5),cat.pos=c(-30,45),cat.dist=c(.04,.04),euler.d=TRUE,filename="E2_DIM_ERa_nearest_gene_overlap_high_res.png",imagetype = "png",height=500,width=500,resolution = 300,lwd=1)
#grid.draw(e2.dim)

combo.e2.dim <- euler(c("E2,ERα" = 4303, "DIM,ERα" = 1349, "E2,ERα&DIM,ERα" = 13289))
tiff("test_qual.tiff",width=800)
plot(combo.e2.dim, counts = TRUE, alpha=0.5,
     fill=c("white", "white"),main="Overlap of genes closest to E2 and DIM ERα binding sites",labels=list(fontsize=5,cex=c(3,3)),quantities=list(cex=c(2,1.5,2)),cat.dist=c(.04,.04))
dev.off()


e2.res=venn.diagram(list(`E2,ERa`=e2.genes,`RES,ERa`=res.genes),filename=NULL,cex=2,fill=c("yellow","green"),main.cex=c(1),cat.cex=c(1.5,1),cat.pos=c(-150,150),cat.dist=c(.04,.025))
grid.draw(e2.res)



combo.e2.res <- euler(c("E2,ERα" = 3245, "RES,ERα" = 1484, "E2,ERα&RES,ERα" = 14347))
plot(combo.e2.res, counts = TRUE, alpha=0.5,
     fill=c("white", "white"),main="Overlap of genes closest to E2 and resveratrol ERα binding sites",labels=list(fontsize=6,cex=4),quantities=list(cex=c(4,4,4)))

tiff("E2_RES_ERa_cloest_genes_overlap.tiff",width=800)
plot(combo.e2.res, counts = TRUE, alpha=0.5,
     fill=c("white", "white"),main="Overlap of genes closest to E2 and RES ERα binding sites",labels=list(fontsize=4,cex=c(3,3)),quantities=list(cex=c(2,1.5,2)),cat.dist=c(.04,.04))
dev.off()



e2.e2.plus.tcdd=venn.diagram(list(`E2,ERa`=e2.genes,`E2+TCDD,ERa`=e2.tcdd.genes),filename=NULL,cex=2,fill=c("yellow","green"),main.cex=c(1),cat.cex=c(1.5,1),cat.pos=c(-150,150),cat.dist=c(.04,.025))
grid.draw(e2.e2.plus.tcdd)


combo.e2.e2.plus.tcdd <- euler(c("E2,ERα" = 1266, "E2+TCDD,ERα" = 2673, "E2,ERα&E2+TCDD,ERα" = 16326))
plot(combo.e2.e2.plus.tcdd, counts = TRUE, font=1, cex=4, alpha=0.5,
     fill=c("white", "white"),main="Overlap of genes closest to ERα ligand induced sites after E2 and E2+TCDD",labels=list(fontsize=15,cex=3),quantities=list(cex=3))




e2.tcdd=venn.diagram(list(`E2,ERa`=e2.genes,`TCDD,ERa`=tcdd.genes),filename=NULL,cex=2,fill=c("yellow","green"),main.cex=c(1),cat.cex=c(1.5,1),cat.pos=c(-150,150),cat.dist=c(.04,.025))
grid.draw(e2.tcdd)


combo.e2.tcdd <- euler(c("E2,ERα" = 17054, "TCDD,ERα" = 70, "E2,ERα&TCDD,ERα" = 538))
plot(combo.e2.tcdd, counts = TRUE, font=1, cex=4, alpha=0.5,
     fill=c("white", "white"),main="Overlap of genes closest to ERα ligand induced sites after E2 and TCDD",labels=list(fontsize=15,cex=3),quantities=list(cex=3))




dim.res=venn.diagram(list(`DIM,ERa`=dim.genes,`RES,ERa`=res.genes),filename=NULL,cex=2,fill=c("yellow","green"),main.cex=c(1),cat.cex=c(1.5,1),cat.pos=c(-150,150),cat.dist=c(.04,.025))
grid.draw(dim.res)



combo.dim.res <- euler(c("DIM,ERα" = 2189, "RES,ERα" = 3382, "DIM,ERα&RES,ERα" = 12449))
plot(combo.dim.res, counts = TRUE, alpha=0.5,
     fill=c("white", "white"),main="Overlap of genes closest to DIM and resveratrol ERα binding sites",labels=list(fontsize=6,cex=4),quantities=list(cex=c(3,3,4)))

tiff("DIM_RES_ERa_cloest_genes_overlap.tiff",width=800)
plot(combo.dim.res, counts = TRUE, alpha=0.5,
     fill=c("white", "white"),main="Overlap of genes closest to DIM and RES ERα binding sites",labels=list(fontsize=4,cex=c(3,3)),quantities=list(cex=c(2,1.5,2)),cat.dist=c(.04,.04))
dev.off()


