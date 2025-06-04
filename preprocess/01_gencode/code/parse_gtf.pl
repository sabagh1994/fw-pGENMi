#!/usr/bin/perl
use List::Util qw(min max);

#chr1    HAVANA  gene    65419   71585   .       +       .       gene_id "ENSG00000186092.5_3"; gene_type "protein_coding"; gene_name "OR4F5"; level 2; havana_gene "OTTHUMG00000001094.3_3"; remap_status "full_contig"; remap_num_mappings 1; remap_target_status "overlap";
my $ifile = $ARGV[0];
my $filter = $ARGV[1];
my $feature = $ARGV[2];

open(IFILE, "<$ifile") or die("Couldn't open $ifile for r");
while(my $line = <IFILE>){
  my @sep = split(/\t/, $line);  
	if ($sep[2] ne $filter) { next;}
	my $chr = $sep[0];
	my $beg = min($sep[3], $sep[4]);
	my $end = max($sep[3], $sep[4]);
	foreach my $field (split(/"; /, $sep[8])){
		if ($field !~ /$feature/) {
			next;
		}
		my ($key, $value) = split(/ /, $field);
		$value =~ s/\"//g;
		print join("\t", ($chr, $beg, $end, $value))."\n";
		last;
	}
}
close(IFILE);
