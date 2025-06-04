#!/bin/sh


## ./xpct.sh ifile q ofile
## Assume the input file is sorted by column $10 ## CHANGE THIS
ifile=$1;
## q is the top/bottom percent of genes to consider by some criteria (column 10)
q=$2;
ofile=$3;
# Get number of genes
n=`wc -l $ifile | cut -d ' ' -f 1`;
# Get q/2. for example, if q=10%, p=5%
p=`perl -se 'print $q/2' -- -q=$q`;
# Get num genes that is p of n. For instance, 5% of 100 for n=100
x=`perl -se 'print int($n*$p)' -- -n=$n -p=$p`;
# Get lines for top x (5) genes with positive 10 column
head -n $x $ifile | awk '{if ($10 > 0) print $_}' > $ofile
# Get lines for bottom x (5) genes with positive 10 column
#tail -n $x $ifile | awk '{if ($10 < 0) print $_}' >> $ofile
