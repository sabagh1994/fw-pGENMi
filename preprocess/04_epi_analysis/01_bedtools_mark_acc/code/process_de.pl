#!/usr/bin/perl
use List::Util qw(min max);
my $ifile = $ARGV[0];
my $argc = $#ARGV + 1;
if ($argc != 1) {
	print "process_de.pl ifile\n";
	exit(1);
}
open(IFILE, "<$ifile") or die("Couldn't open $ifile for reading\n");
while(my $line = <IFILE>){
	chomp($line);
	my @sep = split(/\t/,$line);
	my ($chr1, $beg1, $end1, $name1, $score1, $strand1, $sigval1, $pval1, $qval1, $peak1, $chr2, $beg2, $end2, $name2, $score2, $strand2, $sigval2, $pval2, $qval2, $peak2) = @sep;

	my $begi = max($beg1, $beg2);
	my $endi = min($end1, $end2);
	my $lengthi = $endi - $begi;
	my $namei = join(",", ($name1, $name2));
	my $scorei = $score2 - $score1;
	my $strandi = $strand1;
	my $sigvali = $sigval2 - $sigval1;
	my $pvali = $pval2 - $pval1;
	my $qvali = $qval2 - $qval1;
	my $peaki = $peak2 - $peak1;
	print join("\t",($chr1, $begi, $endi, $namei, $scorei, $strandi, $sigvali, $pvali, $qvali, $peaki))."\n";
}
close(IFILE)
