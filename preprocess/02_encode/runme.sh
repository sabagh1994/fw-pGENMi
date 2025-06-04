#!/bin/sh
#
#
# Goto: encodeproject.org and select Data/Matrix to get to here https://www.encodeproject.org/matrix/?type=Experiment

# OPEN encode_selection.png for correct selection (it still downloads too many files so you will have to parse the metadata)

# Click Download at bottom of page and Download again in the pop-up menu
# Downloaded from encode and renamed: hg19_narrow_broad_HCT-116_TF_encode_files.txt
xargs -n 1 curl -O -L < hg19_narrowbroad_hct116_tfencode.txt  
for i in *.gz; do gunzip $i; done 
# metadata.tsv has info on what TFs each file corresponds to
# Download metadata.tsv into excel
# Column 43: hg19
# Column 3: optimal idr thresholded peaks
# Column 2: bed narrowPeak
# Column 7: HCT-116
# Column 46: released
# Column 17 is the TF with the suffix-human
# Column 1 is the Experiment name - bed suffix
# Column 25 is the date of the experiment
# We will select hg19, optimal idr thresholded peaks, narrowPeak format

awk -F $'\t' '{if ($2 == "bed narrowPeak" && $3 == "optimal idr thresholded peaks" && $7 == "HCT116" && $43 == "hg19" && $46 == "released") { print $1".bed\t"$17"\t"$25}}' metadata.tsv | sed -e "s/-human//g" | sort -u > optimal_idr_files2tf.tsv

# Edit in vim so that ENCFF418WAW.bed CTCF    2017-10-19 is selected among the CTCF experiments (it is the newest)
# Remove third column
# By the end, your file should look like the following:
# ENCFF088WYS.bed JUND
# ENCFF144BSH.bed ZFX
# ENCFF418WAW.bed CTCF
# ENCFF617QEN.bed POLR2AphosphoS5
# ENCFF769NCM.bed EZH2
# ENCFF947SUP.bed EZH2phosphoT487

# create directory called raw and mv all bed files into it
mkdir raw && mv *.bed raw/
mkdir tf_optimal_idr_hg19_HCT-116
while read line
	do
		file=`echo $line | cut -d ' ' -f 1`;
		tf=`echo $line | cut -d ' ' -f 2`;
		cp raw/$file tf_optimal_idr_hg19_HCT-116/$tf;
	done < optimal_idr_files2tf.tsv

# Merge all optimal_idr files into .bed10 file with 4th column replaced by tf name
# You should do this within tf_optimal_idr_hg19_HCT-116 directory
cd tf_optimal_idr_hg19_HCT-116
for tf in *; do
	awk -F $'\t' -v tf=$tf 'BEGIN{OFS=FS} {$4=tf; print $_}' $tf >> ../tfs.bed10
done
cd -

# Sort bed10 file
# Format of the sorted bed file with 10 fields
# chr1    91035   91451   CTCF    612     .       19.55880        -1.00000        5.12429 208
sort -k1,1 -k2,2n tfs.bed10 -o tfs.bed10.sorted
