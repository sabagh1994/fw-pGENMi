#!/bin/sh

run1=1;
run2=1;
run3=1;
run4=1;
run5=1;

perc=2.0
mark="K4me1"

codedir="./mycode";
basedir="./$mark";
final_dir="./final/narrowPeak";

rm -rf $basedir
mkdir -p $basedir
mkdir -p $final_dir

r1dir="../../data/$mark/Rep1";
r2dir="../../data/$mark/Rep2";
p0r1="$r1dir/SW480-A-$mark-5.FCC7C5GACXX_L6_R1_IACCTCA.PE_macs2_peaks.encodePeak";
p0r2="$r2dir/2SW480-A-$mark-5.FCC7WFRACXX_L3_R1_IACCTCA.PE_macs2_peaks.encodePeak";
p6r1="$r1dir/P6-D-$mark-8.FCC7C5GACXX_L7_R1_IGTCGTA.PE_macs2_peaks.encodePeak";
p6r2="$r2dir/2P6-D-$mark-8.FCC7WFRACXX_L4_R1_IGTCGTA.PE_macs2_peaks.encodePeak";


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
cp $basedir/final/de_p0_p6.final.ranked.peak.narrowPeak $final_dir/all/$mark.narrowPeak


# the foloowing steps are for evidence type: TFBS overlapping with mark presence at any stage
# sorting p0/p6 final files based on the peak height => this is later used for 
# evidence type in which the presence of the mark in TFBS matters only, not the 
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

