#!/bin/bash

##################################################################################
# For each model (pgenmi, fwpgenmi), provided by the user                        #
# Step 1 : Perform cross validation on H1 and H0 variants                        #
# Step 2: a) Aggregate cross validation results to find the best                 #
#            hyperparameters (regulatory distance, L2 regularization coeff)      #
#         b) Train on the best hyperparameters                                   #
#         c) Rank the transcription factors                                      #
##################################################################################

modeltype=$1 # user-provided argument \in {pgenmi, fwpgenmi}
SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
PROJPATH=$SCRIPT_PATH/..

source $PROJPATH/activate venv

if [[ $modeltype == pgenmi ]]
then
    # pgenmi
    # Step 1:
    echo "Performing CV on $modeltype"
    python $PROJPATH/fwpgenmi/cv_pgenmi.py --config_path $PROJPATH/configs/01_cfg_cv_pgenmi_H1.yml # H1 variant
    python $PROJPATH/fwpgenmi/cv_pgenmi.py --config_path $PROJPATH/configs/01_cfg_cv_pgenmi_H0.yml # H0 variant
    
    # Step 2:
    echo "Aggregating CV results, training the final model, and ranking TFs for $modeltype"
    python $PROJPATH/fwpgenmi/cvaggr_tr_evidrank.py --config_path $PROJPATH/configs/03_cfg_cvaggr_tr_evidrank_pgenmi.yml
elif [[ $modeltype == fwpgenmi ]]
then
    # fwpgenmi
    # Step 1:
    echo "Performing CV on $modeltype"
    python $PROJPATH/fwpgenmi/cv_fwpgenmi.py --config_path $PROJPATH/configs/02_cfg_cv_fwpgenmi_H1.yml # H1 variant
    python $PROJPATH/fwpgenmi/cv_fwpgenmi.py --config_path $PROJPATH/configs/02_cfg_cv_fwpgenmi_H0.yml # H0 variant

    # Step 2:
    echo "Aggregating CV results, training the final model, and ranking TFs for $modeltype"
    python $PROJPATH/fwpgenmi/cvaggr_tr_evidrank.py --config_path $PROJPATH/configs/04_cfg_cvaggr_tr_evidrank_fwpgenmi.yml
else
    echo "$modeltype is invalid."
fi