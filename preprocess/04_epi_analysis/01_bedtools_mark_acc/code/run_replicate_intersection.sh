#!/bin/sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# rep1 rep2 px path/to/px
mkdir -p $4;
bedtools intersect -a $1 -b $2 -wa -wb > $4/$3.long.narrowPeak
$SCRIPT_DIR/process_intersection.pl $4/$3.long.narrowPeak "$3r1" "$3r2" > $4/$3_intersections.narrowPeak;
sort -k1,1 -k2,2n $4/$3_intersections.narrowPeak -o $4/$3_intersections.narrowPeak;
