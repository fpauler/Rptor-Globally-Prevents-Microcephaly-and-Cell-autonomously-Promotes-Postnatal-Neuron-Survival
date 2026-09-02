the idea is to calculate the coverage of the deleted region for QC
do with bedtools as this allows to deal with split reads!

find * -name "*.Aligned.sortedByCoord.out.bam" | perl del_count/count_reads.pl

this command puts a bed file in each analysis cluster with the counts of reads in this region

