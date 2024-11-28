### This is a representation of the pipeline I used for processing RNA sequencing samples
 
java -jar /Users/siddhaduio.no/Desktop/All_omics_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.fastq.gz /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.fastq.gz.trimmed ILLUMINACLIP:/Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/adapters/TruSeq/TruSeq2-SE.fa:2:30:7 MINLEN:15

/Users/siddhaduio.no/Desktop/All_omics_tools/hisat2-2.2.1/hisat2 -p 8 --dta -x /Users/siddhaduio.no/Desktop/All_omics_tools/hisat2_indexes/hg38_tran/genome_tran -U /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.fastq.gz -S /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.sam --summary-file /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.alignment.summary.txt

/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools view -S -b /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.sam > /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.bam

/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools sort --output-fmt BAM -o /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.sorted.bam /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.bam  

mv /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.sorted.bam /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.bam

/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools index /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.bam

### getting count using htseq. Running union since most approproriate and recommended
htseq-count -f bam -q -m union -s reverse -t exon -i gene_id /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.bam /Users/siddhaduio.no/Desktop/All_omics_tools/Homo_sapiens.GRCh38.104.with.chr.word.gtf > /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/2-3MCF-DMSO_S10_R1_001.htseq.counts.txt



