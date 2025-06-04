#!/bin/bash
set -e

run1=0
if [ $run1 -eq 1 ]
then
	ln -s -f ../01_gencode/protein_coding_genes.list .
fi

tag=$1 # regulatory distance ["10Kb", "50Kb", "200Kb", "1Mb"]
evidtype=$2 # type of the mark evidence being "diffmark": loss/gain of mark/acc, "presmark": presence of mark/acc peak

for i in `cat marks.list`; do mkdir -p binary_feature_$evidtype/$tag/$i; done


eviddir="../04_epi_analysis/02_intersect_tfbs/intersections_$evidtype/$tag/by_mark_$tag/binary/"
datadir="../data/deseq2/"

cd $eviddir
for i in *; do for j in $i/*; do ./code/make_input.pl protein_coding_genes.list $j $datadir/p0vsp6_DESeq_processed > ./binary_feature_$evidtype/$tag/$j; done done

