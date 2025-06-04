#!/bin/bash

#####################################################
###      Making the input for pgenmi/fwpgenmi     ###
#####################################################
#  1. Adding the differential expression p-value    #
#     of genes to the regulatory evidence to        #
#     create the input.	                            #
#  2. The input is separated by TFs                 #
#  3. Combine features of all TFs per evidence type #
#     Generate up vs down-analysis input files      #
#####################################################

rm -rf binary_feature_diffmark
rm -rf binary_feature_presmark

# Creating the input using differental mark peaks (histone mark or accessibility)
./inputgen.sh 10Kb diffmark
./inputgen.sh 50Kb diffmark
./inputgen.sh 200Kb diffmark
./inputgen.sh 1Mb diffmark

# Creating the input using the presence of mark peaks (histone mark or accessibility)
./inputgen.sh 10Kb presmark
./inputgen.sh 50Kb presmark
./inputgen.sh 200Kb presmark
./inputgen.sh 1Mb presmark


# Generate TFBS-only, TFBS_DiffMark, TFBS_DiffAcc, etc features
# aggregate the features for all TFs per evidence type
python ./code/evidence_aggregator.py
