#!/bin/sh

################################################################
### Process the histone mark peaks and accessibility profile ###
################################################################

# Step 1: Intersect the replicates for each mark per stage p0/p6
# ----------------------------------------------------------------------------
# Step 2: Calculate differential peaks (peak_p6 - peak_p0)
# ----------------------------------------------------------------------------
# The regions that are present only in one stage are used as differential peak, i.e., loss/gain of peaks 
# The differential mark at the overlapping regions of p0 and p6 is not considered.
# ----------------------------------------------------------------------------
# Step 3: Generate present peak files for each mark (later used for TFBS overlapping with mark presence)
./K4me1.sh
./K4me3.sh
./K27ac.sh
./K27me3.sh
./ACC.sh

######################################################
### convert narrowpeak to bed format for all marks ###
######################################################
python ./code/narrow2bedconv.py
