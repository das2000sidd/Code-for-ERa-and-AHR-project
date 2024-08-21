#!/bin/bash

## This is a wrapper script that loops over the single end RNA sequencing fastq files and generated individual scripts to process the raw RNA sequencing reads


for fastq_1 in `ls *.fq.gz`
do
			
			
			basename="${fastq_1%.fq.gz}"
			echo 'Fastq is:' $fastq_1
			echo 'Basename:' $basename
			
 
			echo "java -jar /Users/siddhaduio.no/Desktop/All_omics_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 ./$fastq_1 ./$basename.trimmed ILLUMINACLIP:/Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/adapters/TruSeq/TruSeq2-SE.fa:2:30:7 MINLEN:15 \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/hisat2-2.2.1/hisat2 -p 8 --dta -x /Users/siddhaduio.no/Desktop/All_omics_tools/hisat2_indexes/hg38_tran/genome_tran -U ./$basename.fastq.gz -S ./$basename.sam --summary-file ./$basename.alignment.summary.txt \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools view -S -b ./$basename.sam > /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/RNA_Seq/fastq/$basename.bam \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools sort --output-fmt BAM -o ./$basename.sorted.bam ./$basename.bam \n" >> $basename.sh

			echo "mv ./$basename.sorted.bam .q/$basename.bam \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools index ./$basename.bam \n" >> $basename.sh


done
