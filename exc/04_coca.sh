#!/bin/bash

##########################################################################
# Perform COCA and survival analysis                                     #
# Step 1: Generate the input for cluster of cluster analsysis (COCA)     #
#         and survival analysis using different gene sets                #
# Step 2: Fetch the network metadata and get the networks available for  #
#         the provided species ids (9606 for homo sapiens)               #
# Step 3: Perform COCA                                                   #
# Step 4: Aggregate the results of COCA and survival analysis            #
##########################################################################

SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
PROJPATH=$SCRIPT_PATH/..
source $PROJPATH/activate venv

# Step 1:
python $PROJPATH/fwpgenmi/coca_inputgen.py --config_path $PROJPATH/configs/06_cfg_coca_inputgen.yml

# Step 2: 
python $PROJPATH/fwpgenmi/KnowEng.py --NETWORK_DIR_PATH $PROJPATH/results/coca_results/network --species_id 9606

# Step 3: 
python $PROJPATH/fwpgenmi/coca_analysis.py --config_path $PROJPATH/configs/07_cfg_coca_analysis.yml

# Step 4:
python $PROJPATH/fwpgenmi/coca_stats_aggr.py --config_path $PROJPATH/configs/08_cfg_coca_aggr.yml