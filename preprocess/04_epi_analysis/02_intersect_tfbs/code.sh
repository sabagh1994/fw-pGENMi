#!/bin/bash
set -e

run1=1;

# Create soft links to the TFBS files with varying distances once
if [ $run1 -eq 1 ]
then
	addr='../../03_overlap/encode_gencode'
	ln -s -f $addr/tfs_protein_coding_10Kb.bed10 ./
	ln -s -f $addr/tfs_protein_coding_50Kb.bed10 ./
	ln -s -f $addr/tfs_protein_coding_200Kb.bed10 ./
	ln -s -f $addr/tfs_protein_coding_1Mb.bed10 ./
fi

base_dir="../../04_epi_analysis/"

tag=$1 # the regulatory distance, ["1Mb", "200Kb", "50Kb", "10Kb"]
evidtype=$2 # type of the mark evidence being "diffmark": loss/gain of mark/acc, "presmark": presence of mark/acc peak
eviddir=$3 # for diffmark: $base_dir/bedtools_mark_acc/final/bed5/all, presmark: $base_dir/bedtools_mark_acc/Union 
intersections="intersections_"$evidtype

# create intersetions folder and put get_max_by_tf.pl in it 
mkdir -p $intersections/$tag/coords_$tag
mkdir -p $intersections/$tag/gene_sets_$tag
mkdir -p $intersections/$tag/by_mark_$tag

cd $eviddir
for i in *.bed; 
do 
   bedtools intersect -a $base_dir/02_intersect_tfbs/tfs_protein_coding_$tag.bed10 -b $i -wa -wb > $base_dir/02_intersect_tfbs/$intersections/$tag/coords_$tag/$i;
done

## Separate bed files by (TF, Target, Peak Signal) for each mark
cd $base_dir/02_intersect_tfbs/$intersections/$tag/coords_$tag
for i in *.bed; do awk -F $'\t' '{x=gsub(/\|/,"\t",$4); print $4"\t"$14}' $i > $base_dir/02_intersect_tfbs/$intersections/$tag/gene_sets_$tag/$i; done

## Get maximum signals 
cd $base_dir/02_intersect_tfbs/$intersections/$tag/gene_sets_$tag
mkdir -p max_by_tf/continuous max_by_tf/binary
for i in *.bed; do $base_dir/02_intersect_tfbs/get_max_by_tf.pl $i > ./max_by_tf/continuous/$i; done
cd max_by_tf/continuous


if [[ $evidtype -eq "diffmark" ]]
then
  echo "binarizing in diffmark mode"
  ## Convert to pgenmi bit vector with two fields: up down
  for i in *.bed; do awk '{if ($3 < 0) { print $1"\t"$2"\t0\t1"} else { print $1"\t"$2"\t1\t0"} }' $i >> ../binary/$i; done

  ## Make by_mark genesets
  cd ../../
  for i in *.bed; do x=`echo $i | sed -e "s/\.bed//g"`; mkdir -p ../by_mark_$tag/binary/$x; done

  ## Copy to by_mark binary, this separates out each TF to a file >$1
  cd max_by_tf/binary
  for i in *; do x=`echo $i | sed -e "s/\.bed//g"`; awk -v odir=../../../by_mark_$tag/binary/$x/ '{print $2"\t"$3"\t"$4 > odir"/"$1}' $i; done
elif [[ $evidtype -eq "presmark" ]]
then
  echo "binarizing in diffmark mode"
  
  ## Convert to pgenmi bit vector (one element)
  for i in *.bed; do awk '{ print $1"\t"$2"\t1" }' $i >> ../binary/$i; done

  ## Make by_mark genesets
  cd ../../
  for i in *.bed; do x=`echo $i | sed -e "s/\.bed//g"`; mkdir -p ../by_mark_$tag/binary/$x; done

  ## Copy to by_mark binary, this separates out each TF to a file >$1
  cd max_by_tf/binary
  for i in *; do x=`echo $i | sed -e "s/\.bed//g"`; awk -v odir=../../../by_mark_$tag/binary/$x/ '{print $2"\t"$3 > odir"/"$1}' $i; done
else
  echo "invalid evidtype"
fi
