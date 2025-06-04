''' 
    Create gene sets (gene signatures) based on different criterion, which will be later used in the 
    COCA input generator for the COCA analysis.

    Criterion 1: Having a fitted model, get the mediator genes based on the posterior odds ratio, POR 
                (i.e. genes mediating the effect of all TFs).
    
    Criterion 2: Genes mediating the effect of highly ranked TFs (Ratio of post odds ratio, RPOR)
                 Having a fitted model, get the ratio of POR of genes being mediator before and after TF evidence removal 
                 The routine in this code is similar to the one used for getting the evidence significance.

    Criterion 3: DE genes: Use the raw data with information on DE genes (this option is not supported in this script)
'''

import os
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict

from io_cfg import results_dir, input_dir
from io_cfg import make_abspath, read_yaml
from mediator_genes import mediator_gene_run


# -------------------------------- #
#        Get the Arguments         #
# -------------------------------- #

parser = argparse.ArgumentParser()
parser.add_argument('--config_path', type=str, required=True) # e.g., 05_cfg_mediatorgene.yml
args = parser.parse_args()

cfg_path = args.config_path
config_dict = read_yaml.load(cfg_path)
    
# -------------------------------- #
#        Getting the configs       #
# -------------------------------- #

dist = config_dict.get('dist', '50Kb')
dircs = config_dict.get('dircs', ['up', 'down'])
evid_type = config_dict.get('evid_type', 'TFBS_DiffMark')

evid_dir = config_dict.get('evid_dir', 'input')
outdir = config_dict.get('outdir', 'mediator_genes')
evid_sig_outdir = config_dict.get('evid_sig_outdir', 'evid_significance')
trained_model_dir = config_dict.get('trained_model_dir', 'pgenmi_final_results')

top_TF_counts = config_dict.get('top_TF_counts', [5]) # To get the union top X TFs for both analyses, 
                                                      # per gene product of RPORs for the TFs in the union will be computed
top_gene_counts = config_dict.get('top_gene_counts', [10, 20, 50, 70, 100, 200]) # To get the top X genes based on POR or multiplied RPOR criterion

# DE_indir = config_dict.get('DE_indir', '../../preprocess/data/deseq2/p0vsp6_DESeq_processed.csv')

# check for the absolute vs relative path
# ----------------------------------------
evid_dir = make_abspath(evid_dir, input_dir)
outdir = make_abspath(outdir, results_dir)
evid_sig_outdir = make_abspath(evid_sig_outdir, results_dir)
trained_model_dir = make_abspath(trained_model_dir, results_dir)

# -------------------------------------------------------------------------- #
#   Calculate Post Odds Ratio (POR) and Ratio of Post Odds Ratio (RPOR)      #
#          before and after removing TF regulatory evidence                  #
# -------------------------------------------------------------------------- #

for dirc in dircs:
    mediator_gene_run(evid_dir= evid_dir, trained_model_dir=trained_model_dir, outdir= outdir, 
                      output_file=None, evid_type=evid_type, dist=dist, dirc=dirc)
    

# ------------------------------------------------------------------------------------- #
#                         Take the top Genes based on POR                               #
# ------------------------------------------------------------------------------------- #
# 1. Take the top genes based on POR for each direction separately                      #
# 2. Take the union of top genes, ranked by POR, for both directions                    #
#    Justification: top genes in up analysis are ranked at the bottom for down analysis #
#                   and vice versa. To combine both analyses we need to take the union. #
# ------------------------------------------------------------------------------------- #

save_allEvid_genes= False
sep_dic_genes= defaultdict(list) # to save the gene names for up/down
for top in top_gene_counts:
    for dirc in dircs:
        addr= f'{outdir}/{evid_type}/postOddsRatio_org_{dirc}'
        df= pd.read_csv(addr, sep= '\t')
        df.sort_values(by= ['post_OddsRatio_org'], ascending= False, inplace= True)
        
        df_top= df.iloc[:top, :]
        sep_dic_genes[top]= sep_dic_genes[top] + df_top['genes'].tolist()
        save_dir= f'{outdir}/postOddsRatio/Separate/{dirc}'
        os.makedirs(save_dir, exist_ok = True)
        # saving along the post odds ratio
        #df_top.to_csv(f'{save_dir}/postOddsRatio_top{top}_full', sep= '\t', index= False)
        
        # saving the gene names only
        df_top[['genes']].to_csv(f'{save_dir}/postOddsRatio_top{top}', sep= '\t', 
                                 index= False, header= False)
        # save all the genes used for input of pgenmi (protein coding genes + DE genes)
        if not save_allEvid_genes:
            save_allEvid_genes= True
            os.makedirs(f'{outdir}/AllGenes', exist_ok= True)
            df[['genes']].to_csv(f'{outdir}/AllGenes/AllEvidGenes', sep= '\t', header= False, index= False)
            
# taking the union of genes (union of up and down direction) using the post odds ratio
for top in sep_dic_genes:
    lst= list(set(sep_dic_genes[top]))
    save_dir= f'{outdir}/postOddsRatio/Union/'
    os.makedirs(save_dir, exist_ok= True)
    file= open(f'{save_dir}/postOddsRatio_top{top}', 'w')
    file.write('\n'.join(lst))
    file.close()
    

# ------------------------------------------------------------------------------------- #
#                         Take the top Genes based on RPOR                              #
# ------------------------------------------------------------------------------------- #
# 1. Get the ranking of TFs for up and down-analysis then take the union of the top X   #
#    TFs for both direction.                                                            #
# 2. For each gene, multiply the RPOR of TFs in the union set derived from (1)          #
# 3. Rank the genes based on the product of RPORs                                       #
# 4. Save the top x genes based on product of RPOR of top y TFs for each direction.     #
# 5. Save the union of (4) for both directions                                          #
# ------------------------------------------------------------------------------------- #

# getting the candidate TFs
top_TF_dic= {} # holding the candidate top TFs for different number of top TFs
for top in top_TF_counts:
    candid_tfs= []
    for dirc in dircs:
        ranking_dir= f'{evid_sig_outdir}/{evid_type}/'
        #df= pd.read_csv(f'{ranking_dir}/{dirc}/pgenmi_ranking_TFLevel', sep= '\t')
        df= pd.read_csv(f'{ranking_dir}/{dirc}/ranking_TFLevel', sep= '\t')
        df.columns= ['TF', 'loglike', 'delta_llr']
        df.sort_values(by= ['delta_llr'], ascending= False, inplace= True)
        candid_tfs= candid_tfs + df['TF'].tolist()[:top]
    candid_tfs= list(set(candid_tfs))
    top_TF_dic[top]= candid_tfs
    
    # save the top candidate TFs used for mediator genes
    save_dir= f'{outdir}/RatioOfPostOddsRatio/TF_candidates'
    os.makedirs(save_dir, exist_ok =True)
    file= open(f'{save_dir}/topTF_{top}', 'w')
    file.write('\n'.join(candid_tfs))
    file.close()
    
for top_tf in top_TF_counts:
    candid_tfs= top_TF_dic[top]
    for top_gene in top_gene_counts:
        all_genes= []
        for dirc in dircs:
            addr= f'{outdir}/{evid_type}/RatioOfPostOddsRatio_afterTobefore_{dirc}'
            df= pd.read_csv(addr, sep= '\t')
            df['prod']= df[candid_tfs].product(axis=1)
            df.sort_values(by= ['prod'], ascending= True, inplace= True)
            df_top= df.iloc[:top_gene, :]
            genes= df_top.iloc[:,0].tolist()
            all_genes= all_genes + genes
            
            # get the top genes for each direction and save it to a separate file
            save_dir= f'{outdir}/RatioOfPostOddsRatio/Separate/{dirc}'
            os.makedirs(save_dir, exist_ok =True)
            file= open(f'{save_dir}/topGene{top_gene}_topTF{top_tf}', 'w')
            file.write('\n'.join(genes))
            file.close()
            
            
            # also merge them later (do the same as above for post odds ratio)
        save_dir= f'{outdir}/RatioOfPostOddsRatio/Union'
        os.makedirs(save_dir, exist_ok =True)
        all_genes= list(set(all_genes))
         # for each pair of (top tf, top gene) we should save this
        file= open(f'{save_dir}/topGene{top_gene}_topTF{top_tf}', 'w')
        file.write('\n'.join(all_genes))
        file.close()
        
        

# Uncomment to support DE geneset generation:
# -------------------------------------------
# # Criterion 3: DE Genes:
# 
# de_df= pd.read_csv(DE_indir)
# de_df= de_df[['hgnc_symbol', 'pvalue']]

# # dropping the genes with nan names won't have any effect because it 
# # has been handled in preprocessing but I have repeated it here as well :)
# de_df_reduc= de_df.dropna(subset=['hgnc_symbol']) 
# values = {'hgnc_symbol': 'nogene', 'pvalue': 1.}
# # filling the nan p-vapues with 1
# de_df_reduc.fillna(value=values, inplace= True)

# de_df_reduc.sort_values(by= ['pvalue'], inplace= True)

# outdir= 'GeneSets/DEGenes'
# os.makedirs(outdir, exist_ok= True)
# for top in top_gene_counts:
#     genes= de_df_reduc['hgnc_symbol'].tolist()[:top]
#     file= open(f'{outdir}/DEGenes_top{top}', 'w')
#     file.write('\n'.join(genes))
#     file.close()