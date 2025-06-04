#!/usr/bin/perl

my $xfile = $ARGV[0];
my $ifile = $ARGV[1];
my $tfile = $ARGV[2];

my %keep = ();
open(IFILE, "<$xfile") or die("Couldn't open $xfile for r");
while(my $line = <IFILE>){
	chomp($line);
	$keep{$line} = 1;
}
close(IFILE);

open(IFILE, "<$ifile") or die("Couldn't open $ifile for r");
my $length = undef;
while(my $line = <IFILE>){
	chomp($line);
	my @sep = split(/\t/, $line);
	$length = $#sep;
	my $gene = $sep[0];
	if (!defined($keep{$gene})) { next;}
	my $vector = join("\t",@sep[1..$#sep]);
	$geneset{$gene} = $vector;
}
close(IFILE);

# Name    baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
open(TFILE, "<$tfile") or die("Couldn't open $tfile for r");
<TFILE>;
print "PVAL\tA0";
my @n = ();
for (my $i = 1; $i <= $length; $i++){
	push @n, "0";
	print "\tA$i";	
}
my $nstr = join("\t", @n);
print "\n";
while(my $line = <TFILE>){
	chomp($line);
	my ($name, $basemean, $log2foldchange, $lfcse, $pval, $padj, $namerep) = split(/\t/, $line);
	if ($pval eq "NA") {
		$pval = 1;
	}
	if (!defined($keep{$name})) { next;}
	print $name."\t".$pval."\t1\t";
	if (defined($geneset{$name})){
		print $geneset{$name}."\n";
	} else {
		print $nstr."\n";
	}
}
close(TFILE);
