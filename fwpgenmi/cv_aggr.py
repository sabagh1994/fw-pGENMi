import os
import math
import operator
import pandas as pd
from collections import defaultdict
from io_cfg import results_dir, make_abspath

def cv_aggregator(indir='cv_results', outdir='cv_aggr_results', config_dict= None):
    
    '''
    Description: To aggregate the results of cross validation and choose the hyperparameters
                 (regulatory distance and regularization coefficient) to train the final models.
    
    :param indir: (str) path to the cross validation results to aggregate
    :param outdir: (str) path to save the stats of aggregation, which help decide 
                   the best regularization coefficient and regulatory distance for 
                   final model trainings per evidence type.
    :param config_dict: (dict) dictionary containing the list of evidence types 
                        to report the aggregated stats for and the count cross validation repeates.
                        e.g., config_dict = {'evidences': ['TFBS_DiffMark'], 'folds':100}
                        Note that the folds value should match the folds value used for cross validation.
                        
    The following files will be generated and saved:
    1. 'cv_summ_per_cfg': Combined CV results in one file. 
        for all setting tuples (evidence type, H*,distance,direction,regularization coeff, fold#)
        the train and test loglikelihoods are reported.
    2. 'tunning_{direction}': saved in '{outdir}/{evidence_type}', for each evidence type, top 10
       hyperparameter configs (distance, regularization coeff), based on average test llr over the folds, are listed.
       test_llr = H1_llr - H0_llr on the test set.
    3. 'best_config_allEvid_bestUp', 'best_config_allEvid_bestDown': best hyperparam config
       (distance, regularization coeff) for all evidence types are stored separately for each direction of analysis.
    4. 'best_config_allEvid_maxTestLLRSum': For all evidence types, best hyper parameter config 
        is listed based on the sum of averaged test llr for both directions of analysis (up- and down-analyese).
    '''
    
    # check for the absolute vs relative path
    # ----------------------------------------
    indir = make_abspath(indir, results_dir)
    outdir = make_abspath(outdir, results_dir)
    
    evidences = config_dict.get('evidences', ['TFBS_DiffMark'])
    folds = config_dict.get('folds', 100)
    os.makedirs(outdir, exist_ok= True)
    
    #####################################################
    ### Combine the CV results for all evidence types ###
    #####################################################
    
    statfields = ['htag', 'dist', 'dirc', 'regul', 'delta', 'FC', 
                  'train_loglike_total', 'train_loglike', 'test_loglike']
    cv_summ_dict = defaultdict(list)
    for evid_type in evidences:
        # combine the cv runs for the given evidence types
        resdir = f'{indir}/{evid_type}'
        if not(os.path.exists(resdir)):
            print(f'{resdir} does not exist! Continue')
            continue
        for subdir, dirs, files in os.walk(resdir):
            for filename in files:
                if filename.startswith("cfg"):
                    file= open(os.path.join(subdir,filename), 'r')
                    line= file.readlines()[0]
                    stats = line.strip().split('\t')
                    assert len(stats) == len(statfields)
                    for stat_id, stat_val in enumerate(stats):
                        kstat = statfields[stat_id]
                        cv_summ_dict[kstat].append(stat_val)
                    cv_summ_dict['evid_type'].append(evid_type)
        
        column_order = ['evid_type'] + statfields
        cv_summ_df = pd.DataFrame.from_dict(cv_summ_dict)
        cv_summ_df = cv_summ_df[column_order]
        
    # convert stats columns to numeric type
    num_cols = ['train_loglike_total', 'train_loglike', 'test_loglike']
    cv_summ_df[num_cols] = cv_summ_df[num_cols].apply(pd.to_numeric)
    cv_summ_df.to_csv(f'{outdir}/cv_summ_per_cfg', sep='\t', index=False)
    
    ##################################################################
    #           per (evid_type, htag, dirc, regul, delta)            #
    #  compute the sum of train_loglike, test_loglike over the folds #
    ##################################################################
    
    grp_vars = ['evid_type', 'htag', 'dist', 'dirc', 'regul', 'delta']
    target_vars = ['train_loglike', 'test_loglike']
    grouped = cv_summ_df.groupby(grp_vars)
    statsumm_lst = []
    for (name, grp) in grouped:
        assert grp.shape[0] == folds
        evid_type, htag, dist, dirc, regul, delta = name
        tvars = []
        for tvar in target_vars:
            tvar_sum = grp[tvar].sum()
            tvars.append(tvar_sum)
        statsumm_lst.append([evid_type, htag, dist, dirc, regul, delta] + tvars)
        
    statsumm_df = pd.DataFrame(statsumm_lst)
    statsumm_df.columns = grp_vars + target_vars
    
    #########################################################
    #          per (evid_type, dirc, regul, delta)          #
    #  get the average llr, llr = loglike(H1) - loglike(H0) #
    #########################################################
    
    dict_allevid = {}
    grp_vars = ['evid_type', 'dist', 'dirc', 'regul', 'delta']
    target_vars = ['train_loglike', 'test_loglike']
    grouped = statsumm_df.groupby(grp_vars)

    for (name, grp) in grouped:
        if (grp.shape[0] == 1) and (grp['htag'].squeeze()=='H0'):
            print(f'Corressponding H1 run does not exist for the H0 setting {name}! Exiting!')
            break
        assert (grp.shape[0] == 1) or (grp.shape[0] == 2)
        evid_type, dist, dirc, regul, delta = name

        if not evid_type in dict_allevid.keys():
            dict_allevid[evid_type] = defaultdict(list)

        H1_test = grp[grp['htag']=='H1']['test_loglike'].squeeze()
        H1_train = grp[grp['htag']=='H1']['train_loglike'].squeeze()

        # condition on all grp_vars except regul, regul is set to 0 for H0 runs
        df_h0 = statsumm_df[(statsumm_df['htag']=='H0') & (statsumm_df['evid_type']==evid_type) &
                            (statsumm_df['dist']==dist) & (statsumm_df['dirc']==dirc) & 
                            (statsumm_df['delta']==delta)]
        assert df_h0.shape[0] == 1, 'more than one regularization coefficient has been tried for H0'
        H0_test = df_h0['test_loglike'].squeeze()
        H0_train = df_h0['train_loglike'].squeeze()

        test_llr_ave = (H1_test - H0_test)/folds
        train_llr_ave = (H1_train - H0_train)/folds

        dist_no = 1000 if 'Mb' in dist else int(dist.split('Kb')[0])
        dict_allevid[evid_type][dirc].append([dist_no, regul, test_llr_ave, train_llr_ave, dist])
        

    final_configs_top10= {}
    for evid_type in dict_allevid.keys():
        final_configs_top10[evid_type]={}
        for dirc in dict_allevid[evid_type].keys(): #dircs:
            lst = dict_allevid[evid_type][dirc]
            # take the top 10 high average test llr configs
            lst_sorted_test_llr = sorted(lst, key=operator.itemgetter(2), reverse=True)[0:10]
            final_configs_top10[evid_type][dirc] = lst_sorted_test_llr

    #final_configs_top10
    for evid_type in evidences:
        for dirc in final_configs_top10[evid_type].keys(): #dircs:
            nested_lst = final_configs_top10[evid_type][dirc]
            df = pd.DataFrame(nested_lst, columns= ['distance', 'lambda', 
                                                    'ave_test_llr', 'ave_train_llr', 'dist_str'])
            path = f'{outdir}/{evid_type}/'
            os.makedirs(path, exist_ok= True)
            df.to_csv(f'{path}/tunning_{dirc}')



    # for each evidence, take the best config for up and use it for down later 
    # for each evidence, take the best config for down and use it for up later 
    file_up = open(f'{outdir}/best_config_allEvid_bestUp', 'w')
    file_down = open(f'{outdir}/best_config_allEvid_bestDown', 'w')
    for evid_type in evidences:
        if 'up' in final_configs_top10[evid_type].keys():
            [dist_no, regul, test_llr_ave, train_llr_ave, dist] = final_configs_top10[evid_type]['up'][0]
            file_up.write(f'{evid_type}\t{dist}\t{regul}\t{test_llr_ave}\t{train_llr_ave}\n')

        if 'down' in final_configs_top10[evid_type].keys():
            [dist_no, regul, test_llr_ave, train_llr_ave, dist] = final_configs_top10[evid_type]['down'][0]
            file_down.write(f'{evid_type}\t{dist}\t{regul}\t{test_llr_ave}\t{train_llr_ave}\n')

    file_up.close()
    file_down.close()
    
    
    # for each evidence type sum average test/train LLR for both directions
    # then take the config that has the highest overall test LLR 
    final_configs_all_summed= {}
    for evid_type in dict_allevid.keys():
        lst_up = dict_allevid[evid_type]['up'] # list of lists
        lst_down = dict_allevid[evid_type]['down'] # list of lists 

        lst_up = sorted(lst_up, key=operator.itemgetter(0), reverse=True)
        lst_up = sorted(lst_up, key=operator.itemgetter(1), reverse=True)

        lst_down = sorted(lst_down, key=operator.itemgetter(0), reverse=True)
        lst_down = sorted(lst_down, key=operator.itemgetter(1), reverse=True)

        summed_lst = []
        for i in range(len(lst_up)):
            [dist_no_up, regul_up, test_llr_ave_up, train_llr_ave_up, dist_up] = lst_up[i]
            [dist_no_down, regul_down, test_llr_ave_down, train_llr_ave_down, dist_down] = lst_down[i]

            if (dist_no_up != dist_no_down) or (regul_up != regul_down):
                print('distance or regul mismatch!')
                break

            test_llr_ave_sum = test_llr_ave_up + test_llr_ave_down
            train_llr_ave_sum = train_llr_ave_up + train_llr_ave_down

            summed_lst.append([dist_no_up, regul_up, test_llr_ave_sum, train_llr_ave_sum, dist_up])
        # sort based on summed average test llr for both directions
        summed_lst = sorted(summed_lst, key=operator.itemgetter(2), reverse=True)
        # take the highest value
        final_configs_all_summed[evid_type] = summed_lst[0]
        
        # saving all the summed stats (summed average test/train llr) for the evidence type => no filtering
        df_summed_all = pd.DataFrame(summed_lst, columns=["distance_no", "lambda", "test_llr_ave_sum", 
                                                          "train_llr_ave_sum", "evid_distance"])
        os.makedirs(f'{outdir}/all_summed_stat/', exist_ok=True)
        df_summed_all.to_csv(f'{outdir}/all_summed_stat/{evid_type}_allSummed', sep= '\t')
        
    file_summ = open(f'{outdir}/best_config_allEvid_maxTestLLRSum', 'w')
    for evid_type in evidences:
        [dist_no, regul, test_llr_ave, train_llr_ave, dist] = final_configs_all_summed[evid_type]
        file_summ.write(f'{evid_type}\t{dist}\t{regul}\t{test_llr_ave}\t{train_llr_ave}\n')
    file_summ.close()