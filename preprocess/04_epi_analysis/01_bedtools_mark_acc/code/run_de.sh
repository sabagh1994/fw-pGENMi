#!/bin/sh

mkdir -p $4
name1="p0";
name2="p6";

#bedtools intersect -a $1 -b $2 -wa -wb > $4/$3.long.narrowPeak
bedtools intersect -a $1 -b $2 -v > $4/$name1.exclusive.narrowPeak
bedtools intersect -a $2 -b $1 -v > $4/$name2.exclusive.narrowPeak
#$basedir/mycode/process_de.pl $4/$3.long.narrowPeak > $4/$3_intersections.narrowPeak;
#sort -k1,1 -k2,2n $4/$3_intersections.narrowPeak -o $4/$3_intersections.narrowPeak;
sort -k1,1 -k2,2n $4/$name1.exclusive.narrowPeak -o $4/$name1.exclusive.narrowPeak;
sort -k1,1 -k2,2n $4/$name2.exclusive.narrowPeak -o $4/$name2.exclusive.narrowPeak;
#sort -k1,1 -k2,2n $4/$3.long.narrowPeak -o $4/$3.long.narrowPeak;
awk '{print $1,$2,$3,$4,-1*$5,$6,-1*$7,-1*$8,-1*$9,-1*$10}' $4/$name1.exclusive.narrowPeak | tr ' ' '\t' > blah && mv blah $4/$name1.exclusive.narrowPeak;

cat $4/$name1.exclusive.narrowPeak > $4/$3.final.narrowPeak
cat $4/$name2.exclusive.narrowPeak >> $4/$3.final.narrowPeak
#cat $4/$3_intersections.narrowPeak >> $4/$3.final.narrowPeak
sort -k1,1 -k2,2n $4/$3.final.narrowPeak -o $4/$3.final.narrowPeak;
sort -grk 10 $4/$3.final.narrowPeak -o $4/$3.final.ranked.peak.narrowPeak;
