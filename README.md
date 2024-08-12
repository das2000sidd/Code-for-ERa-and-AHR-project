# Code-for-ERa-and-AHR-project
This repository holds the prototype of script used for processing raw sequencing data and subsequent downstream analysis for generating fiugres used in the paper looking at crosstalk between arylhydrocarbon receptor and estrogen receptor alpha project.
The GEO repository for the raw data for this project is deposited in GSE232235.

The script ChIPSeq_Sample_processing.sh holds the code used for processing raw single end ChIP sequencing experiment FASTQ files which included trimming Illumina TruSeq adapter using Trimmomatic (v.0.39) followed by alignment to the reference genome using bowtie2 to generate SAM files. The SAM files were converted to bam files using samtools view flag followed sorting of the BAM file by samtools. Subsequently the BAM files were indexed using samtools index flag.

The script RNASeq_Sample_processing.sh holds the code used for processing raw single end RNA sequencing experiment FASTQ files which included trimming Illumina TruSeq adapter using Trimmomatic (v.0.39) followed by alignment to the reference genome using bowtie2 to generate SAM files. The SAM files were converted to bam files using samtools view flag followed sorting of the BAM file by samtools. Subsequently the BAM files were indexed using samtools index flag. htseq2-count was used to generate raw count matrices for the RNA seq data.

The script differential_expression_analysis.R has the script that uses the count file generated from RNASeq_Sample_processing.sh for differential expression analysis.

The script get_distance_to_TSS_ERa_peaks.R takes as input a set of peaks and calculates and generates distance to TSS plot for those set of peaks using the GenomicDistributions and annoatates the peaks using the ChIPseeker package. 

The script UpSet_plot_overlapping_DEGs.R generates an upset plot of overlap among DEG's of various RNA sequencing differential expression analysis of relevance. 



The published paper is at: https://www.mdpi.com/1422-0067/24/19/14578
