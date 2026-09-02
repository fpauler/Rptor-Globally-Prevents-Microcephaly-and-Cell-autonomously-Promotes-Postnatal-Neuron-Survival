# analyses read counts in genomic intervals using bedtools

#input bam files comes from pipe - usually with find * -name "*.Aligned.sortedByCoord.out.bam"

while (<STDIN>) {
	chomp($_);

	@parts = split(/\//, $_);
	$folder = "$parts[0]";
	for ($i=1; $i<$#parts; $i++) {
		$folder = $folder."/".$parts[$i];
	}

	$cmd = "bedtools intersect -c -b $_ -a Rptor_exon6.bed -split -bed > $folder/Rptor_exon6.counts.out";
	print "$cmd\n";
	system($cmd);
}



