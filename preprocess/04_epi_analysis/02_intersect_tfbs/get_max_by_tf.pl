#!/usr/bin/perl

my $ifile = $ARGV[0];
my %tf2gene=();
open(IFILE,"<$ifile") or die("Die r $ifile\n");
while(my $line = <IFILE>){
	chomp($line);
	my ($tf, $gene, $score) = split(/\s+/,$line);
	if (!defined($tf2gene{$tf}{$gene}) || abs($tf2gene{$tf}{$gene}) < $score){
		$tf2gene{$tf}{$gene} = $score;
	}
}
close(IFILE);

foreach my $tf (keys %tf2gene){
	foreach my $gene (keys %{$tf2gene{$tf}}){
		my $score = $tf2gene{$tf}{$gene};
		print "$tf\t$gene\t$score\n";
	}
}
