import torch
import random
import torch.nn as nn
import numpy as np
import os

import factpgenmi_torch
from factpgenmi_torch import get_file_dims, load_data, PGM_LL_modif, estep, mstep, partition_data


def factpgenmi_torch_final_run(indir= '../../input', outdir= 'factpgenmi_final_results', regul_= 10, dirc= 'up',
                           evid_type= 'TFBS_DiffMark', dist= '50Kb', perc= 0.8, delta= 0.001, is_H0= False,
                           mode= 'all', tf_dim= 20, mark_dim= 8, output_file= None):

    assert mode in ['train', 'all'], 'mode is not valid: %r' % mode
    #assert evid_type in ['TFBS_DiffACC', 'TFBS_DiffMark',
    #                     'TFBS_DiffMarkAggr', 'TFBS_only', 'TFBS_PresMark'], 'evidence type is not valid: %r' % evid_type

    if not is_H0:
        regul= regul_
    else:
        regul= 0

    if isinstance(output_file, str):
        output_file_obj = open(output_file, 'w')
    else:
        output_file_obj = output_file


    if is_H0:
        out_path= f'{outdir}/{evid_type}/{dirc}/H0'
        input_file= f'{indir}/{evid_type}/{dist}/{dirc}/H0_{dirc}'
    else:
        out_path= f'{outdir}/{evid_type}/{dirc}/H1'
        input_file= f'{indir}/{evid_type}/{dist}/{dirc}/H1_{dirc}'

    os.makedirs(out_path, exist_ok=True)

    num_genes, num_columns = get_file_dims(input_file, output_file=output_file_obj)
    num_features = num_columns - 2;
    genes, data, evidence = load_data(input_file, num_features)

    if isinstance(output_file, str):
        output_file_obj.close()

    if is_H0:
        evidence= np.reshape(evidence, (num_genes, 1))

    # converting numpy to tensor
    if mode == 'train':
        (train_inds, test_inds)= partition_data(num_genes, perc= perc)
        data_train= data[train_inds]
        evidence_train= evidence[train_inds, :]
        data_test= data[test_inds]
        evidence_test= evidence[test_inds, :]

        evidence_tensor_test= torch.from_numpy(evidence_test).double()
        data_tensor_test= torch.from_numpy(data_test).double()

        evidence_tensor= torch.from_numpy(evidence_train).double()
        data_tensor= torch.from_numpy(data_train).double()
    elif mode == 'all':
        data_tensor= torch.from_numpy(data).double()
        evidence_tensor= torch.from_numpy(evidence).double()

    # this part is extra:
    if len(evidence_tensor.shape) ==1:
        evidence_tensor= torch.reshape(evidence_tensor, (num_genes, 1))
    num_features= evidence_tensor.shape[1]

    # model creation:
    model= PGM_LL_modif(tf_dim= tf_dim, mark_dim= mark_dim)
    # parameter initializations:
    for name, p in model.named_parameters():
        if name == 'pgm_w_tf' or name == 'pgm_w_m' or name == 'pgm_bias':
            torch.nn.init.xavier_uniform_(p) #_ does it inplace
        elif name == 'alpha':
            torch.nn.init.uniform_(p)
        else:
            print(f'unknown parameter! {name}')

    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    min_steps= 10
    max_steps= 15000
    last_n= 100
    delta_loglike= delta+1
    i= 1
    while i<= min_steps or (delta_loglike > delta and i<=max_steps):
        loglike_org, loglike, new_a, post_zg1, pgm_w= estep(model, evidence_tensor,
                                                            data_tensor, regul, is_H0)
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


    torch.save(model.state_dict(), f'{out_path}/saved_model')
    if mode == 'train': # needs to be adjasted (it is not used here so it is not adjasted)
        with torch.no_grad():
            test_loglike, _, _, _, _= estep(model, evidence_tensor_test,
                                            data_tensor_test, regul, is_H0)

        file= open(f'{args.addr}/{regul}_{delta}_{FC}_{args.rank}', 'w')
        file.write(f'{regul}\t{delta}\t{FC}\t{loglike_new}\t{loglike_org}\t{test_loglike}\t{perc}\n')
        file.close()
    else:
        loglike_np= loglike_org.detach().numpy().reshape(-1)
        pgm_w_np= np.array(pgm_w.detach()).reshape(-1)
        pgm_alpha= np.array(model.alpha.detach()).reshape(-1)
        params_merged= np.concatenate((pgm_w_np, pgm_alpha), axis= 0)
        params_merged= np.concatenate((loglike_np, params_merged), axis= 0)

        pgm_w_tf_np= np.array(model.pgm_w_tf.detach()).reshape(-1)
        pgm_w_m_np= np.array(model.pgm_w_m.detach()).reshape(-1)
        tf_mark_weights= np.concatenate((pgm_w_tf_np, pgm_w_m_np), axis= 0)

        np.savetxt(f'{out_path}/factpgenmi_final_stat', params_merged, delimiter= '\t')
        np.savetxt(f'{out_path}/factpgenmi_TFMarkWeights', tf_mark_weights, delimiter= '\t')
        return loglike_np
