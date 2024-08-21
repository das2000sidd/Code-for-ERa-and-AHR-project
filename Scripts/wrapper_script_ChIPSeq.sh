#!/bin/bash

## This is a wrapper script that loops over the ChIP seq fastq files and generated individual scripts to process the raw ChIP sequencing reads

for fastq_1 in `ls *.fastq.gz`
do
			
			
			basename="${fastq_1%.fastq.gz}"
			echo 'Fastq is:' $fastq_1
			echo 'Basename:' $basename
			
 
			echo "java -jar /Users/siddhaduio.no/Desktop/All_omics_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 ./$fastq_1 ./$basename.clean.fq.trimmed ILLUMINACLIP:/Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/adapters/TruSeq/TruSeq2-SE.fa:2:30:7 MINLEN:15 \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/bowtie2-2.4.4-macos-x86_64/bowtie2 -p 64 -q -x /Users/siddhaduio.no/Desktop/All_omics_tools/bowtie2-2.4.4-macos-x86_64/human_hg38_index -U ./$basename.clean.fq.trimmed -S ./$basename.clean.sam 2> /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/$basename.bowtie.log \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools view -S -b ./$basename.clean.sam > ./$basename.clean.bam \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools sort -o ./$basename.sorted.clean.bam ./$basename.clean.bam \n" >> $basename.sh

			echo "mv ./$basename.sorted.clean.bam ./$basename.clean.bam \n" >> $basename.sh

			echo "/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools index ./$basename.clean.bam \n" >> $basename.sh


done
