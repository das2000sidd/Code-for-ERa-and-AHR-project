### This is a representation of the pipeline I used for processing ChIP sequencing samples

 
java -jar /Users/siddhaduio.no/Desktop/All_omics_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.fq.gz /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.fq.trimmed ILLUMINACLIP:/Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/adapters/TruSeq/TruSeq2-SE.fa:2:30:7 MINLEN:15

/Users/siddhaduio.no/Desktop/All_omics_tools/bowtie2-2.4.4-macos-x86_64/bowtie2 -p 64 -q -x /Users/siddhaduio.no/Desktop/All_omics_tools/bowtie2-2.4.4-macos-x86_64/human_hg38_index -U /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.fq.trimmed -S /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.sam 2> /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.bowtie.log

/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools view -S -b /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.sam > /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.bam

/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools sort -o /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.sorted.clean.bam /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.bam

mv /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.sorted.clean.bam /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.bam

/Users/siddhaduio.no/Desktop/All_omics_tools/samtools-1.14/samtools index /Users/siddhaduio.no/Desktop/PhD_Project_related/Data_from_Jason/ChIP_Seq/raw_reads/DIM_AHR_Rep_1.clean.bam


