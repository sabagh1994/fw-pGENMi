#!/usr/bin/perl
use List::Util qw(max min);

my $ifile = $ARGV[0];

# chr1    9907    10642   DPF2,NBN,ZBTB40,ZBTB48,ZNF146,ZNF239    968.5   .       118.5501325     -1      3.1827325       99.05555555     chr1    0       65419   OR4F5
open(IFILE, "<$ifile") or die("Couldn't open $ifile for r");
while(my $line = <IFILE>){
  chomp($line);
  my ($chr, $p1b, $p1e, $name1, $s5, undef, $s7, $s8, $s9, undef, $p2b, $p2e, $name2) = split(/\s+/, $line);
  my $start = max($p1b, $p2b);
  my $stop = min($p1e, $p2e);
  my $name = $name1."|".$name2;
  if (!defined($name2)){
    $name = $name1;
  }
  print "$chr\t$start\t$stop\t$name\t$s5\t.\t$s7\t$s8\t$s9\n";
}
close(IFILE);
