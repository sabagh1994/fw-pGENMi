import os
import torch
import numpy as np
import pandas as pd

from io_cfg import results_dir, input_dir, make_abspath
from fwpgenmi import PGM_LL, PGM_LL_FW, estep
from fwpgenmi import get_file_dims, load_data

def evidence_ranker(evid_dir='input', trained_model_dir='fwpgenmi_final_results', outdir='evid_significance',
                    output_file= None, evid_type='TFBS_DiffMark', dist='50Kb', dirc='up', 
                    tf_dim=None, mark_dim=None, model_type='fwpgenmi'):
    
    '''
    Description: Rank the TFs based on their contribution to the fit of model to the data
    
    :param evid_dir: (str) rootpath to the input evidence files, used for training
    :param evid_type: (str) type of the evidence e.g., "TFBS_DiffMark", "TFBS_only", etc. 
    :param dist: (str) regulatory distance \in {"10Kb", "50Kb", "200Kb", "1Mb"}
    :param dirc: (str) direction of the analysis, up or down
    :param model_type: (str) model_type \in {"pgenmi", "fwpgenmi"}
    :param trained_model_dir: (str) path to the saved model parameters
                              e.g., trained_model_dir \in {'fwpgenmi_final_results', 'pgenmi_final_results'}
    :param outdir: (str) path to save the results of TF rankings and posterior probs
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
    # where the ranking is saved and the posterior probabilities
    out_path= f'{outdir}/{evid_type}/{dirc}'
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

    # loadind in the trained model parameters
    saved_model_path = f'{trained_model_dir}/{evid_type}/{dirc}/H1/saved_model'
    model.load_state_dict(torch.load(saved_model_path))

    data_tensor= torch.from_numpy(data).double()
    evidence_tensor= torch.from_numpy(evidence).double()
    num_features= evidence_tensor.shape[1]
    
    # because we only care about the original loglike it does not matter 
    # what regularization is used
    loglik_org, _, _, post_zg1 = estep(model, evidence_tensor, data_tensor, regul=0, is_H0=False)
    
    df_post_zg1= pd.DataFrame(post_zg1.detach().numpy(), index= genes)
    #df_post_zg1.to_csv(f'{out_path}/{model_type}_post_zg1', sep= '\t')
    df_post_zg1.to_csv(f'{out_path}/post_zg1', sep= '\t')
    df_data= pd.read_csv(input_file, sep= '\t')
    feat_labels= list(df_data.columns[1:])
    
    
    TF_names= []
    if evid_type != 'TFBS_only':
        for i, evid in enumerate(feat_labels):
            if i==0:
                continue
            TF_names.append(evid.split('_')[0])
        TF_names= list(set(TF_names))
    else:
        TF_names= feat_labels[1:]

    tf_ind_dict= {}
    for tf in TF_names:
        #tf_inds= [index+1 for index, value in enumerate(feat_labels[1:]) if value.split('_')[0] == tf]
        tf_inds = [index for index, value in enumerate(feat_labels[0:]) if value.split('_')[0] == tf]
        tf_ind_dict[tf] = tf_inds
        
    dic= {}
    for tf in tf_ind_dict.keys():
        tf_inds= tf_ind_dict[tf]
        evidence_tensor_cp= evidence_tensor.clone()
        evidence_tensor_cp[:, tf_inds]= 0
        loglik_after, _, _, _ = estep(model, evidence_tensor_cp, data_tensor, regul=0, is_H0= False)
        delta_log= loglik_org - loglik_after
        dic[tf]= [np.asscalar(loglik_after.detach().numpy()), np.asscalar(delta_log.detach().numpy())]
       
    df = pd.DataFrame.from_dict(dic).transpose()
    df= df.sort_values(1, ascending= False)
    #df.to_csv(f'{out_path}/{model_type}_ranking_TFLevel', sep= '\t')
    df.to_csv(f'{out_path}/ranking_TFLevel', sep= '\t')
