#!/bin/sh

# The codes and comments come from the anlaysis of differential histone mark peaks 
# so some of the comments should be adjusted

run1=1;
run2=1;
run3=1;
run4=1;
run5=1;

perc=2.0
mark="ACC"

codedir="./code";
basedir="./$mark";
final_dir="./final/narrowPeak";

rm -rf $basedir
mkdir -p $basedir
mkdir -p $final_dir

ACCdir="../../data/ATAcSeq/"
p0r1="$ACCdir/M0_1101/M0_1101.FCHCTNYBBXX_L4_ITAAGGCGA__hg19_q30_macs_peaks.narrowPeak";
p0r2="$ACCdir/M0_1102/M0_1102.FCHCTNYBBXX_L4_IAGGCAGAA__hg19_q30_macs_peaks.narrowPeak";
p6r1="$ACCdir/M6_1101/M6_1101.FCHCTNYBBXX_L4_ICGTACTAG__hg19_q30_macs_peaks.narrowPeak";
p6r2="$ACCdir/M6_1102/M6_1102.FCHCTNYBBXX_L4_ITCCTGAGC__hg19_q30_macs_peaks.narrowPeak";

## Step 1: Merge Replicates: Produce x/x_intersections.narrowPeak files
if [ $run1 -eq 1 ]
then
        $codedir/run_replicate_intersection.sh $p0r1 $p0r2 p0 $basedir/p0
        $codedir/run_replicate_intersection.sh $p6r1 $p6r2 p6 $basedir/p6
fi


if [ $run2 -eq 1 ]
then
        $codedir/aggregate_merged.pl $basedir/p0/p0_intersections.narrowPeak > $basedir/p0/p0_intersection.final;
        $codedir/aggregate_merged.pl $basedir/p6/p6_intersections.narrowPeak > $basedir/p6/p6_intersection.final;
fi

if [ $run3 -eq 1 ]
then
        $codedir/run_de.sh $basedir/p0/p0_intersection.final $basedir/p6/p6_intersection.final de_p0_p6 $basedir/final;
fi

# copy the results of de peak analysis
mkdir -p $final_dir/all
cp $basedir/final/de_p0_p6.final.ranked.peak.narrowPeak $final_dir/all/DiffAcc.narrowPeak


# the foloowing steps are for evidence type: TFBS overlapping with mark presence at any stage

# sorting p0/p6 final files based on the peak height => this is later used for 
# evidence type in which the presence of the ACC in TFBS matters only, not the 
# presence of a differential peak
if [ $run4 -eq 1 ]
then
	sort -rgk 10 $basedir/p0/p0_intersection.final -o $basedir/p0/p0_intersection.final.peaksorted
	sort -rgk 10 $basedir/p6/p6_intersection.final -o $basedir/p6/p6_intersection.final.peaksorted
	# take the top/bottom m% of mark sites based on the 10th column, i.e. peak heightv
	$codedir/xpct.sh $basedir/p6/p6_intersection.final.peaksorted $perc $basedir/p6/p6_intersection.final.peaksorted.pct
	$codedir/xpct.sh $basedir/p0/p0_intersection.final.peaksorted $perc $basedir/p0/p0_intersection.final.peaksorted.pct
fi
# merge p0 and p6 final files then sort them and take top perc, e.g., 10% of the sorted merged file
# this step is being done to generate evidence types that do not count for diff peaks and
# only consider the presence of mark peak within TFBS in either of the stages p0/p6
if [ $run5 -eq 1 ]
then
	mkdir -p $basedir/union
	cp $basedir/p0/p0_intersection.final $basedir/union/p0_p6.union.NarrowPeak
	cat $basedir/p6/p6_intersection.final >> $basedir/union/p0_p6.union.NarrowPeak
	sort -rgk 10 $basedir/union/p0_p6.union.NarrowPeak -o $basedir/union/p0_p6.union.NarrowPeak.sorted
	$codedir/xpct.sh $basedir/union/p0_p6.union.NarrowPeak.sorted $perc $basedir/union/p0_p6.union.NarrowPeak.sorted.pct
	mkdir -p ./Union
	awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $10}' $basedir/union/p0_p6.union.NarrowPeak.sorted.pct > ./Union/"$mark.bed"
fi


