## This is a prototype of a script I used for ChIP seq peak calling

export PATH=$PATH:/Users/siddhaduio.no/Desktop/All_omics_tools/homer/bin/


findMotifsGenome.pl E2_ER_top_1000_peaks.bed hg38 E2_ERa_top_1000_peaks_homer_motifs
findMotifsGenome.pl DIM_ER_top_1000_peaks.bed hg38 DIM_ERa_top_1000_peaks_homer_motifs
findMotifsGenome.pl RES_ER_top_1000_peaks.bed hg38 RES_ERa_top_1000_peaks_homer_motifs

