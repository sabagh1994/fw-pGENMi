#!/usr/bin/perl
use List::Util qw(min max);
my $ifile = $ARGV[0];
my $prefix1 = $ARGV[1];
my $prefix2 = $ARGV[2];
my $argc = $#ARGV + 1;
if ($argc != 3) {
	print "process.pl ifile prefix1 prefix2\n";
	exit(1);
}
open(IFILE, "<$ifile") or die("Couldn't open $ifile for reading\n");
while(my $line = <IFILE>){
	chomp($line);
	my @sep = split(/\t/,$line);
	if ($prefix ne "-" || $prefix2 ne "-") {
		$sep[3] = $prefix1."_".$sep[3];
		$sep[13] = $prefix2."_".$sep[13];
	}
 	my $num_fields = $#sep + 1;
	if ($num_fields != 20) {
		next;
	}
	my ($chr1, $beg1, $end1, $name1, $score1, $strand1, $sigval1, $pval1, $qval1, $peak1, $chr2, $beg2, $end2, $name2, $score2, $strand2, $sigval2, $pval2, $qval2, $peak2) = @sep;

	my $begi = max($beg1, $beg2);
	my $endi = min($end1, $end2);
	my $lengthi = $endi - $begi;
	my $namei = join(",", ($name1, $name2));
	my $scorei = join(",", ($score1, $score2));
	my $strandi = join(",", ($strand1, $strand2));
	my $sigvali = join(",", ($sigval1, $sigval2));
	my $pvali = join(",", ($pval1, $pval2));
	my $qvali = join(",", ($qval1, $qval2));
	my $peaki = join(",", ($peak1, $peak2));
	print join("\t",($chr1, $begi, $endi, $namei, $scorei, $strandi, $sigvali, $pvali, $qvali, $peaki))."\n";
	#my $begl = min($beg1, $beg2);
	#if ($lengthi != $overlap) {
	#	print "ERROR: $lengthi != $overlap $begi $endi\n";
	#	exit(1);
	#}
	#my $endr = max($end1, $end2);
	#if ($begi != $begl) {
	#	my $i = $beg1 == $begl ? 3 : 13;
	#	my $j = $beg1 == $begl ? 9 : 19;
	#	print "chr1\t$begl\t$begi\t";
	#	print join("\t",@sep[$i..$j])."\n";
	#}
	#if ($endi != $endr) {
	#	$i = $end1 == $endr ? 3 : 13;
	#	$j = $end1 == $endr ? 9 : 19;
	#	print "chr1\t$endi\t$endr\t";
	#	print join("\t",@sep[$i..$j])."\n";
	#}
}
close(IFILE)
