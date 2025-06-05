import os
import torch
import random
import numpy as np
import pandas as pd

from fwpgenmi import PGM_LL, PGM_LL_FW, estep
from fwpgenmi import get_file_dims, load_data
from io_cfg import results_dir, input_dir, make_abspath

def posteriorOddsRatio(post_zg1):
    '''
    Compute the posterior odds ratio p/(1-p)
    
    :param post_zg1: (torch.tensor) posterior probability 
                     of genes being the mediator (Pr(Zg=1))
    :returns: (numpy.array) posterior odds ratio for all genes Pr(Zg=1)/Pr(Zg=0)
    '''
    # the input has tensor type
    post_zg1_np= post_zg1.detach().numpy()
    post_zg0_np= 1 - post_zg1_np

    zero_inds= np.where(post_zg0_np == 0)[0]
    post_zg0_np[zero_inds] = 10**-50

    post_zg1_np= post_zg1_np.reshape(-1,1)
    post_zg0_np= post_zg0_np.reshape(-1,1)
    post_odds_ratio= post_zg1_np/post_zg0_np
    
    return post_odds_ratio


def mediator_gene_run(evid_dir='input', trained_model_dir='pgenmi_final_results', outdir='mediator_genes', 
                      output_file=None, evid_type='TFBS_DiffMark', dist='200Kb', dirc='up', 
                      tf_dim=None, mark_dim=None, model_type='pgenmi'):
    '''
    Compute the posterior odds ratio (POR) for all genes and ratio of post odds ratio (RPOR)
    for all (TF, gene) pairs, for each direction of analysis.
    
    :param evid_dir: (str) rootpath to the evidence input files.
    :param outdir: (str) path to the store the results:
    :param trained_model_dir: (str) path to the saved pretrained model parameters.
    :param model_type: (str) model_type \in {"pgenmi", "fwpgenmi"}
    :param evid_type: (str) evidence type used to access the evidence input file as 
                       well as constructing the path to save the results.
                       evid_type \in {"TFBS_DiffMark", "TFBS_only", "TFBS_DiffMarkAggr", etc}
    :param dist: (str) regulatory distance, used to access the evidence input file
                  dist \in {"10Kb", "50Kb", "200Kb", "1Mb"}
    :param dirc: (str) direction of analyis \in {"up", "down"}
    :param tf_dim: (int) default to None. Number of TFs. Required for fwpgenmi.
    :param mark_dim: (int) default to None. Number of marks*dim of mark representation. Required for fwpgenmi.
    '''
    
    # check for the absolute vs relative path
    # ----------------------------------------
    evid_dir = make_abspath(evid_dir, input_dir)
    outdir = make_abspath(outdir, results_dir)
    trained_model_dir = make_abspath(trained_model_dir, results_dir)
    # ----------------------------------------

    # evidence file
    input_file= f'{evid_dir}/{evid_type}/{dist}/{dirc}/H1_{dirc}'
    # Path to save the rankings and posterior probs
    out_path= f'{outdir}/{evid_type}'
    os.makedirs(out_path, exist_ok= True)

    num_genes, num_columns = get_file_dims(input_file, output_file=output_file)
    num_features = num_columns - 2
    genes, data, evidence = load_data(input_file, num_features);
    
    # Instantiating the model:
    if model_type == 'fwpgenmi':
        assert not(tf_dim is None) and not(mark_dim is None)
        model = PGM_LL_FW(tf_dim=tf_dim, mark_dim=mark_dim)
    elif model_type == 'pgenmi':
        model= PGM_LL(in_dim=num_features, out_dim=1)
    else:
        raise NotImplementedError(f'model type {model_type} is invalid! provide fwpgenmi or pgenmi as model type.')
    
    # Loadig the trained model parameters
    saved_model_path= f'{trained_model_dir}/{evid_type}/{dirc}/H1/saved_model'
    model.load_state_dict(torch.load(saved_model_path)) 

    data_tensor= torch.from_numpy(data).double()
    evidence_tensor= torch.from_numpy(evidence).double() 
    num_features= evidence_tensor.shape[1]

    # We only care about the loglikelihood of data - regularization/prior term does not count
    loglik_org, _, _, post_zg1 = estep(model, evidence_tensor, data_tensor, is_H0=False)
    post_OddsRatio_np= posteriorOddsRatio(post_zg1)

    # save the gene names along with the posterior odds ratio
    df_postOddsRatio_org = pd.DataFrame({'genes': genes.reshape(-1), 'post_OddsRatio_org': post_OddsRatio_np.reshape(-1)})
    df_postOddsRatio_org.sort_values(by=['post_OddsRatio_org'], inplace=True, ascending=False)
    df_postOddsRatio_org.to_csv(f'{out_path}/postOddsRatio_org_{dirc}', sep='\t', index=False)
    
    df_data= pd.read_csv(input_file, sep= '\t')
    feat_labels= list(df_data.columns[1:])


    TF_names= []
    if evid_type != 'TFBS_only':
        for i, evid in enumerate(feat_labels):
            if i==0:
                continue
            TF_names.append(evid.split('_')[0])
        TF_names = list(set(TF_names))
    else:
        TF_names = feat_labels[1:]

    tf_ind_dict= {}
    for tf in TF_names:
        tf_inds = [index for index, value in enumerate(feat_labels[0:]) if value.split('_')[0] == tf]
        tf_ind_dict[tf] = tf_inds


    dic= {} # holding the drop in post odds ratio after TF evidence removal
    for tf in tf_ind_dict.keys():
        tf_inds = tf_ind_dict[tf]
        evidence_tensor_cp = evidence_tensor.clone()
        evidence_tensor_cp[:, tf_inds] = 0

        # Compute posterior probability of genes being mediators
        loglik_after, _, _, post_zg1_after = estep(model, evidence_tensor_cp, data_tensor, regul=0, is_H0= False)
        delta_log = loglik_org - loglik_after

        post_OddsRatio_after_np = posteriorOddsRatio(post_zg1_after)
        # ratio of posterior probs
        ratio_of_postOdds = post_OddsRatio_after_np/post_OddsRatio_np
        dic[tf] = list(ratio_of_postOdds.reshape(-1))

    df = pd.DataFrame.from_dict(dic)
    df.index= genes
    df.to_csv(f'{out_path}/RatioOfPostOddsRatio_afterTobefore_{dirc}', sep= '\t')