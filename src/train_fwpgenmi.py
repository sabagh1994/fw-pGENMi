import os
import torch
import random
import numpy as np
import torch.nn as nn

from fwpgenmi import PGM_LL_FW, estep, mstep
from fwpgenmi import get_file_dims, load_data, partition_data
from io_cfg import results_dir, input_dir, make_abspath


def train_wrapper(indir='input', outdir='fwpgenmi_final_results', regul_=10, dirc='up', 
                  evid_type='TFBS_DiffMark', dist='50Kb', delta= 0.001, is_H0=False, tf_dim=20, 
                  mark_dim=8, output_file=None, save_model=True):
    '''
    Wrapper function to train pgenmi model
    
    :param indir: (str) root path to the evidence files. The input files for each evidence 
                  type (evid_type) should be stored at {indir}/{evid_type}/{dist}/{dirc}/H*_{dirc}
    :param evid_type: (str) evidence type, used to create path to the input file, 
                      e.g., "TFBS_DiffMark", "TFBS_only", "TFBS_PresMark"
    :param dirc: (str) direction of analysis, needed to create the path to the input file
    :param dist: (str) regulatory distance used for the input file, needed to create the path 
                  to the input file, e.g., "10Kb", "50Kb", "200Kb", "1Mb"
    :param is_H0: (bool) whether to include regulatory evidence (False, H1) or not (True, H0)
    :param outdir: (str) root path to saving the results of training
    :param regul_: (float) the regularization coefficient for L2 term in the training objective
    :param tf_dim: (int) number of transcription factors
    
    :param mark_dim: (int) (number of marks) x (the dimension for mark representation), e.g., 
                     Example 1: four histone marks with 2-dimensional one-hot encoding 
                                for the loss/gain in the peak => [mark_dim = 4x2 = 8]
                     Example 2: Accessibility peak with loss/gain of peak, [mark_dim = 1x2 = 2]
                     Example 3: four histone mark with presence of marks in either 
                                stage p0/p6 [mark_dim = 4x1 = 4]

    :param save_model: (bool) save the trained model
    :param output_file: (str or file object) name of the file or a file object used for 
                        the log of reading the evidence input file
                        
    :returns: (numpy arrray) loglikelihood of the fitted model on the training data
    '''
    
    # check for the absolute vs relative path
    # ----------------------------------------
    indir = make_abspath(indir, input_dir)
    outdir = make_abspath(outdir, results_dir)
    # ----------------------------------------

    regul = regul_ if not is_H0 else 0
    if is_H0:
        out_path = f'{outdir}/{evid_type}/{dirc}/H0'
        input_file = f'{indir}/{evid_type}/{dist}/{dirc}/H0_{dirc}' 
    else:
        out_path = f'{outdir}/{evid_type}/{dirc}/H1'
        input_file = f'{indir}/{evid_type}/{dist}/{dirc}/H1_{dirc}'

    os.makedirs(out_path, exist_ok=True)  
    
    # train fwpgenmi model
    cfg_dict = dict(is_H0=is_H0, delta=delta, regul=regul, mark_dim=mark_dim, tf_dim=tf_dim)
    loglike_org, model = train(input_file, output_file, cfg_dict)

    if save_model:
        torch.save(model.state_dict(), f'{out_path}/saved_model')
    
    loglike_np = loglike_org.detach().numpy().reshape(-1)
    pgm_w_np = np.array(model.pgm_w.detach()).reshape(-1)
    pgm_alpha = np.array(model.alpha.detach()).reshape(-1)
    params_merged = np.concatenate((pgm_w_np, pgm_alpha), axis= 0)
    params_merged = np.concatenate((loglike_np, params_merged), axis= 0)

    pgm_w_tf_np = np.array(model.pgm_w_tf.detach()).reshape(-1)
    pgm_w_m_np = np.array(model.pgm_w_m.detach()).reshape(-1)
    tf_mark_weights = np.concatenate((pgm_w_tf_np, pgm_w_m_np), axis= 0)

    np.savetxt(f'{out_path}/fwpgenmi_final_stat', params_merged, delimiter= '\t')
    np.savetxt(f'{out_path}/fwpgenmi_TFMarkWeights', tf_mark_weights, delimiter= '\t')
    return loglike_np


def train(input_file, output_file, cfg_dict):
    '''
    Train routine for the pgenmi model
    
    :param input_file: (str) full path to the input file
    :param output_file: (str or file object) path to the file or a file object used for 
                        the log of reading the evidence input file
    :param cfg_dict: (dict) dictionary containing the parameters used for training
                     and setting up the model 
                     
    :returns: (torch.tensor) loglikelihood of the fitted model on the data
              (torch.nn.Module) trained model
    '''
    
    is_H0 = cfg_dict['is_H0']
    delta = cfg_dict['delta']
    regul = cfg_dict['regul']
    tf_dim = cfg_dict['tf_dim']
    mark_dim = cfg_dict['mark_dim']
    
    output_file_obj = open(output_file, 'w') if isinstance(output_file, str) else output_file
    
    num_genes, num_columns = get_file_dims(input_file, output_file=output_file_obj)
    num_features = num_columns - 2
    genes, data, evidence = load_data(input_file, num_features)
    
    if isinstance(output_file, str):
        output_file_obj.close()
    
    if is_H0:
        evidence= np.reshape(evidence, (num_genes, 1))
    
    # converting numpy to tensor
    data_tensor = torch.from_numpy(data).double()
    evidence_tensor = torch.from_numpy(evidence).double()

    # this part is extra:
    if len(evidence_tensor.shape) ==1:
        evidence_tensor = torch.reshape(evidence_tensor, (num_genes, 1))
    num_features = evidence_tensor.shape[1]

    # -------------------------------- #
    #       Instantiate the model      #
    # -------------------------------- #
    if not is_H0:
        # instantiate fwpgenmi model
        # ---------------------------
        assert num_features == (tf_dim*mark_dim)+1, f'num_features {num_features}!=(tf_dim {tf_dim})*(mark_dim {mark_dim})+(bias 1)'
        # model creation:
        model= PGM_LL_FW(tf_dim=tf_dim, mark_dim=mark_dim)
        # parameter initializations:
        for name, p in model.named_parameters():
            if name == 'pgm_w_tf' or name == 'pgm_w_m' or name == 'pgm_bias':
                torch.nn.init.xavier_uniform_(p)
            elif name == 'alpha':
                torch.nn.init.uniform_(p)
            else:
                print(f'unknown parameter! {name}')
    else:
        # instantiate pgenmi model
        # -------------------------
        # model creation:
        model = PGM_LL(in_dim=num_features, out_dim=1)
        # parameter initializations:
        for name, p in model.named_parameters():
            if name == 'pgm_w':
                torch.nn.init.xavier_uniform_(p)
            elif name == 'alpha':
                torch.nn.init.uniform_(p)
            else:
                print(f'unknown parameter! {name}')
        
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    min_steps = 10
    max_steps = 15000
    last_n = 100
    delta_loglike = delta+1
    i= 1
    while i<= min_steps or (delta_loglike > delta and i<=max_steps):
        loglike_org, loglike, new_a, post_zg1 = estep(model, evidence_tensor, data_tensor, regul, is_H0)
        if i==1:
            loglike_prev = loglike.clone().detach()
            loglike_init = loglike.clone().detach()
        elif (i % last_n == 0):
            loglike_new= loglike.clone().detach()
            with torch.no_grad():
                delta_loglike = loglike_new - loglike_prev
                loglike_prev = loglike_new
        mstep(model, optimizer, loglike, new_a)
        i+=1
        
    return loglike_org, model