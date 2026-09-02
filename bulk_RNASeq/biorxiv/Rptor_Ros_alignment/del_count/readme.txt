the idea is to calculate the coverage of the delted region for QC
do with bedtools as this allows to deal with split reads!

find * -name "*.Aligned.sortedByCoord.out.bam" | perl del_count/count_reads.pl

this command puts a bed file in each analysis cluster with the counts of reads in this region
Rptor ex6 coordinates are based on GRCm39/mm39 
Note: Rptor is too lowly expressed to be detected by this shallow sequencing

