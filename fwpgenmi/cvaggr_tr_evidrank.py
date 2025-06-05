import argparse
import pandas as pd

from cv_aggr import cv_aggregator
from evid_rank import evidence_ranker
import fwpgenmi.train_pgenmi, fwpgenmi.train_fwpgenmi

from fwpgenmi.io_cfg import results_dir, input_dir
from fwpgenmi.io_cfg import make_abspath, read_yaml

# -------------------------------- #
#        Get the Arguments         #
# -------------------------------- #

parser = argparse.ArgumentParser()
parser.add_argument('--config_path', type=str, required=True) # e.g., 03_cfg_cvaggr_tr_evidrank_pgenmi.yml
args = parser.parse_args()

# -------------------------------- #
#        Get the Configs           #
# -------------------------------- #

cfg_path = args.config_path
config_dict = read_yaml.load(cfg_path)

evidences = config_dict.get('evidences', ['TFBS_DiffMark'])
dircs = config_dict.get('dircs', ['up', 'down'])
model_type = config_dict.get('model_type') # "pgenmi", "fwpgenmi"

evid_dir = config_dict.get('evid_dir', 'input')
aggr_outdir = config_dict.get('aggr_outdir', 'cv_aggr_results')
cv_results_dir = config_dict.get('cv_results_dir', 'cv_results')
model_outdir = config_dict.get('model_outdir', 'fwpgenmi_final_results')
evid_sig_outdir = config_dict.get('evid_sig_outdir', 'evid_significance')

tf_dim, mark_dim = None, None
if model_type == "fwpgenmi":
    tf_dim = config_dict.get('tf_dim') # e.g., 20
    mark_dim = config_dict.get('mark_dim') # e.g. 8
    assert not(tf_dim is None) and not(mark_dim is None)

# check for the absolute vs relative path
# ----------------------------------------
evid_dir = make_abspath(evid_dir, input_dir)
aggr_outdir = make_abspath(aggr_outdir, results_dir)
model_outdir = make_abspath(model_outdir, results_dir)
cv_results_dir = make_abspath(cv_results_dir, results_dir)
evid_sig_outdir = make_abspath(evid_sig_outdir, results_dir)

# -------------------------------------------- #
# aggregating the results of cross validation  #
# -------------------------------------------- #
cv_aggregator(indir=cv_results_dir, outdir=aggr_outdir, config_dict=config_dict)

# best_config_allEvid_maxTestLLRSum or best_config_allEvid_bestUp for fwpgenmi
# best_config_allEvid_maxTestLLRSum for pgenmi
input_file = f'{aggr_outdir}/best_config_allEvid_maxTestLLRSum'

# get the best config for each evidence type and train the model on the entire dataset
file = open(input_file, 'r')
final_model_stats = {}

for line in file:
    evid_type, dist, regul, _, _ = line.strip().split('\t')
    final_model_stats[evid_type] = {}
    
    #print(evid_type)
    for dirc in dircs:
        if model_type == "fwpgenmi":
            H1_loglike_org = train_fwpgenmi.train_wrapper(indir=evid_dir, outdir=model_outdir, regul_=float(regul),
                                                          dirc=dirc, evid_type=evid_type, dist=dist, is_H0=False, 
                                                          tf_dim=tf_dim, mark_dim=mark_dim, output_file=None)
        elif model_type == "pgenmi":
            H1_loglike_org = train_pgenmi.train_wrapper(indir=evid_dir, outdir=model_outdir, regul_=float(regul), 
                                                        dirc=dirc, evid_type=evid_type, dist=dist, is_H0=False, output_file=None)
        else:
            raise NotImplementedError(f'invalid model type {model_type}. Supported model types are pgenmi and fwpgenmi')
        
        H0_loglike_org = train_pgenmi.train_wrapper(indir=evid_dir, outdir=model_outdir, regul_=float(regul), 
                                                    dirc=dirc, evid_type=evid_type, dist=dist, is_H0=True, output_file=None)
        LLR = H1_loglike_org - H0_loglike_org
        final_model_stats[evid_type][dirc] = [LLR[0], H1_loglike_org[0], H0_loglike_org[0]]  
        #print(LLR[0], H1_loglike_org[0], H0_loglike_org[0])
file.close()

cols = pd.MultiIndex.from_product([dircs, ['LLR', 'H1_loglike', 'H0_loglike']])
# you need to merge the down and up for each evidence type in the dict
merged_results_dict = {}
for evid_type in evidences:
    merged_results_dict[evid_type] = final_model_stats[evid_type]['up'] + final_model_stats[evid_type]['down']
    
# save the stats for the training on the entire dataset
df = pd.DataFrame(merged_results_dict).transpose()
df.columns = cols
df.to_csv(f'{model_outdir}/final_summary_stats', sep= '\t')

# ------------------------- #
#   Ranking the evidence    #
# ------------------------- #
file = open(input_file, 'r')
for line in file:
    evid_type, dist, _, _, _ = line.strip().split('\t')
    
    for dirc in dircs:
        evidence_ranker(evid_dir=evid_dir, trained_model_dir=model_outdir, outdir=evid_sig_outdir, 
                        output_file=None, evid_type=evid_type, dist=dist, dirc=dirc,
                        tf_dim=tf_dim, mark_dim=mark_dim, model_type=model_type)
file.close()
