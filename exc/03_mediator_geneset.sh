#!/bin/bash

####################################################
#  To get mediator genes by Post Odds Ratio (POR)  # 
#      and Ratio of Post Odds Ratio (RPOR)         #
####################################################

SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
PROJPATH=$SCRIPT_PATH/..

source $PROJPATH/activate venv
python $PROJPATH/src/geneset_gen.py --config_path $PROJPATH/configs/05_cfg_mediatorgene.yml
