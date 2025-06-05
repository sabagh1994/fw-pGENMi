import os
import sys
import argparse

import torch
import random
import numpy as np
import torch.nn as nn
from itertools import product

from fwpgenmi.io_cfg import results_dir, input_dir
from fwpgenmi.io_cfg import make_abspath, read_yaml

from fwpgenmi.fwpgenmi import PGM_LL, estep, mstep
from fwpgenmi.fwpgenmi import get_file_dims, load_data, partition_data

# -------------------------------- #
#        Get the Arguments         #
# -------------------------------- #
parser = argparse.ArgumentParser()
parser.add_argument('--config_path', default='./config_crossValid_pgenmi.json', type=str)
parser.add_argument('--rank', action="store", type=int, default=-1, help="job id in the batch job submission.")
parser.add_argument('--parallel', default=False, action='store_true',
                    help="running configs in parallel. If this option is True then the rank/job_id in the batch job submission should be provided using --rank argument.")
args = parser.parse_args()

cfg_path = args.config_path
config_dict = read_yaml.load(cfg_path)

folds = config_dict.get('folds', 20)
perc = config_dict.get('perc', 0.8)
is_H0 = config_dict.get('is_H0', False)
indir = config_dict.get('indir', 'input')
outdir = config_dict.get('outdir', 'cv_results')
dircs = config_dict.get('dircs', ['up', 'down'])
evid_type = config_dict.get('evid_type', 'TFBS_DiffMark')
distances = config_dict.get('distances', ['10Kb', '50Kb', '200Kb', '1Mb'])

# check if relative or abs path is provided
# ------------------------------------------
indir = make_abspath(indir, input_dir)
outdir = make_abspath(outdir, results_dir)
# ------------------------------------------

output_file= None
if not is_H0:
    lambdas = config_dict.get('lambdas', [100, 10, 1, 0.1, 0.01, 0.001, 0])
else:
    lambdas= [0]

curr_rank= -1
delta= 0.001

cfgs = list(product(distances, lambdas, dircs))

# -------------------------------- #
#       Running the Configs        #
# -------------------------------- #

for cfg_id, cfg in enumerate(cfgs):
    dist, regul, dirc = cfg

    for FC in range(folds):
        if args.parallel:
            curr_rank = curr_rank + 1
            if not (curr_rank == args.rank or args.rank == -1):
                continue
        random.seed((FC+1)*10)

        # generate directories here
        if is_H0:
            #final_outdir= f'{outdir}/{evid_type}/{dist}/{dirc}/H0'
            final_outdir = f'{outdir}/{evid_type}/H0'
            os.makedirs(final_outdir, exist_ok=True)
            input_file = f'{indir}/{evid_type}/{dist}/{dirc}/H0_{dirc}' 
        else:
            #final_outdir= f'{outdir}/{evid_type}/{dist}/{dirc}/H1'
            final_outdir = f'{outdir}/{evid_type}/H1'
            os.makedirs(final_outdir, exist_ok=True)
            input_file = f'{indir}/{evid_type}/{dist}/{dirc}/H1_{dirc}' 


        if isinstance(output_file, str):
            output_file_obj = open(output_file, 'w')
        else:
            output_file_obj = output_file

        num_genes, num_columns = get_file_dims(input_file, output_file=output_file_obj)
        num_features = num_columns - 2
        genes, data, evidence = load_data(input_file, num_features)
        if isinstance(output_file, str):
            output_file_obj.close()

        if is_H0:
            evidence = np.reshape(evidence, (num_genes, 1))
            
        # partition the data
        (train_inds, test_inds) = partition_data(num_genes, perc=perc)
        data_train = data[train_inds]
        evidence_train = evidence[train_inds, :]
        data_test = data[test_inds]
        evidence_test = evidence[test_inds, :]

        # converting numpy to tensor
        evidence_tensor_test = torch.from_numpy(evidence_test).double()
        data_tensor_test = torch.from_numpy(data_test).double()

        evidence_tensor = torch.from_numpy(evidence_train).double()
        data_tensor = torch.from_numpy(data_train).double()

        # extra not needed!
        if len(evidence_tensor.shape) ==1:
            evidence_tensor = torch.reshape(evidence_tensor, (num_genes, 1))
        num_features = evidence_tensor.shape[1]

        # model creation:
        model = PGM_LL(in_dim=num_features, out_dim=1)
        # parameter initializations:
        for name, p in model.named_parameters():
            if name == 'pgm_w':
                torch.nn.init.xavier_uniform_(p)
            elif name == 'alpha':
                torch.nn.init.uniform_(p)
            else:
                print(f'unknown parameter! {name} - not supported')

        optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

        min_steps= 10
        max_steps= 15000
        last_n= 100
        delta_loglike= delta+1
        i= 1
        while i<= min_steps or (delta_loglike > delta and i<=max_steps):
            train_loglike, loglike, new_a, post_zg1 = estep(model, evidence_tensor, data_tensor, regul, is_H0)
            if i==1:
                loglike_prev= loglike.clone().detach()
                loglike_init= loglike.clone().detach()
            elif (i % last_n == 0):
                loglike_new= loglike.clone().detach()
                with torch.no_grad():
                    delta_loglike= loglike_new-loglike_prev
                    loglike_prev= loglike_new
            mstep(model, optimizer, loglike, new_a)
            i+=1

        with torch.no_grad():
            test_loglike, _, _, _ = estep(model, evidence_tensor_test, data_tensor_test, regul, is_H0)

        htag = 'H0' if is_H0 else 'H1'
        file = open(f'{final_outdir}/cfg{cfg_id}', 'w')
        file.write(f'{htag}\t{dist}\t{dirc}\t{regul}\t{delta}\t{FC}\t{loglike_new}\t{train_loglike}\t{test_loglike}\n')
        file.close()