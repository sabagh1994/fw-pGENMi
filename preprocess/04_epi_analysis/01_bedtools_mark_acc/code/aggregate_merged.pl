#!/usr/bin/perl

my $ifile = $ARGV[0];
my $argc = $#ARGV + 1;
if ($argc != 1) {
	print "aggregated_merged.pl ifile\n";
	exit(1);
}
open(IFILE, "<$ifile") or die("Couldn't open $ifile for reading\n");
while(my $line = <IFILE>){
	chomp($line);
	my ($chr, $beg, $end, $names, $scores, $strands, $sigvals, $pvals, $qvals, $peaks) = split(/\t/, $line);
	my $score = average($scores);
	my $strand = ".";
	my $sigval = average($sigvals);
	my $pval = average($pvals);
	my $qval = average($qvals);
	my $peak = average($peaks);
	print join("\t", ($chr, $beg, $end, $names, $score, $strand, $sigval, $pval, $qval, $peak))."\n";
}
close(IFILE);

sub average {
	my $csv = shift;
	my @sep = split(/,/,$csv);
	my $num = $#sep + 1;
	if ($num == 1) {
		return $csv;
	}
	my $sum = $sep[0];
	for (my $i = 1; $i < $num; $i++){
		$sum = $sum + $sep[$i];
	}
	my $mean = $sum/(1.0*$num);
	return $mean;
}
