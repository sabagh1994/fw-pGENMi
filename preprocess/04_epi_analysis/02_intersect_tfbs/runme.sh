#!/bin/sh

###########################################################################################
###           Diff/Pres mark overlapping with TFBS evidence generation:                 ###
###########################################################################################
# Mark stands for either a histone mark or accessibility       			          #
# 1. take the TFBS file with different distance thresholds (10Kb, 50Kb, 200Kb, and 1Mb)   #
# 2. intersect TFBS with diff/pres mark sites						  #
# 3. binarize the evidence: diffmark, gain in the peak => 10, loss in the peak => 01      #
#                           presmark, one regulatory bit is used 0/1                      #
# Results will be stored at intersections_diffmark and intersections_presmark             #
#											  #
###########################################################################################

rm -rf intersections_diffmark
rm -rf intersections_presmark

# Intersect TFBS with differential histone mark/Accessibility peaks 
./code.sh 10Kb diffmark ../../epi_analysis/01_bedtools_mark_acc/final/bed5/all
./code.sh 50Kb diffmark ../../epi_analysis/01_bedtools_mark_acc/final/bed5/all
./code.sh 200Kb diffmark ../../epi_analysis/01_bedtools_mark_acc/final/bed5/all
./code.sh 1Mb diffmark ../../epi_analysis/01_bedtools_mark_acc/final/bed5/all

# Intersect TFBS with present histone mark/Accessibility peaks in either of the stages p0 or p6
./code.sh 10Kb presmark ../../epi_analysis/01_bedtools_mark_acc/Union
./code.sh 50Kb presmark ../../epi_analysis/01_bedtools_mark_acc/Union
./code.sh 200Kb presmark ../../epi_analysis/01_bedtools_mark_acc/Union
./code.sh 1Mb presmark ../../epi_analysis/01_bedtools_mark_acc/Union
