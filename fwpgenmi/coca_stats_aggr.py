import argparse
import pandas as pd
from io_cfg import results_dir, input_dir
from io_cfg import make_abspath, read_yaml

# -------------------------------- #
#        Get the Arguments         #
# -------------------------------- #
parser = argparse.ArgumentParser()
parser.add_argument('--config_path', type=str, required=True) # e.g., 08_cfg_coca_aggr.yml
args = parser.parse_args()

# -------------------------------- #
#         Get the Configs          #
# -------------------------------- #
cfg_path = args.config_path
config_dict = read_yaml.load(cfg_path)
    
gene_counts= config_dict.get('gene_counts', [50, 70, 100])
outdir= config_dict.get('outdir', "COCA_results")
outdir = make_abspath(outdir, results_dir)

if "input_data_dirs" in config_dict.keys():
    input_data_dirs= config_dict['input_data_dirs']
    for i in range(len(input_data_dirs)):
        input_datadir = input_data_dirs[i]
        input_datadir = make_abspath(input_datadir, results_dir)
        input_data_dirs[i] = input_datadir
else:
    #gene_no= [60000, 16000]
    input_data_dirs= [f'{results_dir}/coca_results/coca_input/mediator_genes/AllGenes/AllEvidGenes']
    for i in gene_counts:
        #gene_no= gene_no + 7*[i]
        input_data_dirs.append(f'{results_dir}/coca_results/coca_input/mediator_genes/postOddsRatio/Separate/down/postOddsRatio_top{i}')
        input_data_dirs.append(f'{results_dir}/coca_results/coca_input/mediator_genes/postOddsRatio/Separate/up/postOddsRatio_top{i}')
        input_data_dirs.append(f'{results_dir}/coca_results/coca_input/mediator_genes/RatioOfPostOddsRatio/Separate/down/topGene{i}_topTF5')
        input_data_dirs.append(f'{results_dir}/coca_results/coca_input/mediator_genes/RatioOfPostOddsRatio/Separate/up/topGene{i}_topTF5')

        input_data_dirs.append(f'{results_dir}/coca_results/coca_input/mediator_genes/postOddsRatio/Union/postOddsRatio_top{i}')
        input_data_dirs.append(f'{results_dir}/coca_results/coca_input/mediator_genes/RatioOfPostOddsRatio/Union/topGene{i}_topTF5')
        #input_data_dirs.append(f'coca_results/coca_input/DEGenes/DEGenes_top{i}')
    
aggr_lst= []
for i, subdir in enumerate(input_data_dirs):
    addr = f'{subdir}/results/stats.txt'
    df = pd.read_csv(addr, sep= '\t')
    
    # filtering the dataframe if needed
    #df= df[df.coca_type == 'allOmics']
    #clstr_size= [3, 4]
    #df= df[df.num_clusters.isin(clstr_size)]
    #smoothing= [0.3, 0.8]
    #df= df[df.network_influence.isin(smoothing)]
    #nw_type= ['./network/Property/9606/PPI_complex/9606.PPI_complex.edge', 
    #          './network/Property/9606/enrichr_pathway/9606.enrichr_pathway.edge']
    #df= df[df.network_type.isin(nw_type)]
    
    df.sort_values(by=['survival_analysis_pval'], inplace=True)
    lst= df.iloc[0,:].tolist()
    lst[3]= lst[3].split('/')[-1]
    lst.append('_'.join(subdir.split('/')[2:]))
    #lst.append(gene_no[i])
    aggr_lst.append(lst)
    
df_aggr = pd.DataFrame.from_records(aggr_lst, 
                                    columns= ['pval', 'cluster_no', 'nw_smoothing', 
                                              'nw_type', 'folder_no', 'coca_type', 
                                              'input_type']) #, 'gene_no'])
df_aggr= df_aggr.sort_values(by='pval').reset_index(drop=True)

# removing the Union and postOddsRatio Criterion
# df_aggr= df_aggr[~df_aggr.input_type.str.contains("_postOddsRatio_")]
# df_aggr= df_aggr[~df_aggr.input_type.str.contains("_Union_")]

cols_reorder= ['input_type', 'coca_type', 'folder_no', #'gene_no', 
               'cluster_no', 'nw_smoothing', 'nw_type', 'pval']
df_aggr= df_aggr[cols_reorder]
df_aggr.to_csv(f'{outdir}/aggregated_stat', sep= '\t', index= False)