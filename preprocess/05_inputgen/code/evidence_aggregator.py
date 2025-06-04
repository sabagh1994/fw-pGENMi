import os
import math
import numpy as np
import pandas as pd
import collections

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__)) # preprocess/05_inputgen/code
evid_rootdir = os.path.realpath(f'{SCRIPT_DIR}/../evidences') # preprocess/05_inputgen/evidences
rootdir = os.path.realpath(f'{SCRIPT_DIR}/../..') # preprocess 

DESeq_addr= os.path.join(rootdir, "data/deseq2/p0vsp6_DESeq_processed") #'../data/deseq2/p0vsp6_DESeq_processed'
de_gene_df= pd.read_csv(DESeq_addr, sep= '\t')

# used for TFBS-only evidence generation
pgdir = os.path.join(rootdir, "01_gencode/protein_coding_genes.list")
protein_gene_df= pd.read_csv(pgdir, header= None)
protein_coding_genes= list(protein_gene_df[0])

dists= ['10Kb','200Kb','50Kb','1Mb']

##################################
###      TFBS only evidence    ###
##################################

out_path = f'{evid_rootdir}/TFBS_only/'
os.makedirs(out_path, exist_ok=True)
    
for dist in dists:
    tfbs_gene_inter = os.path.join(rootdir, f"03_overlap/encode_gencode/tfs_protein_coding_{dist}.bed10")
    tf_gene_df = pd.read_csv(tfbs_gene_inter, header= None, sep= '\t')
    # build a dictionary with genes as its keys and the TFs targetting each genes as the values
    gene_dict = collections.OrderedDict()
    tf_dict = collections.OrderedDict()
    tf_id = 0
    # generating a dict of genes and their nearby tfs
    for index, row in tf_gene_df.iterrows():
        pair = row[3].split('|')
        tf = pair[0]
        gene = pair[1]
        if not tf in tf_dict.keys():
            tf_id = tf_id+1
            tf_dict[tf] = tf_id
        if not gene in gene_dict.keys():
            lst = []
            lst.append(tf)
            gene_dict[gene] = lst
        else:
            if not tf in gene_dict[gene]:
                gene_dict[gene].append(tf)
                
    # Note that: 
    # 1. there might be genes with nan p-values you should replace those p-values with 1
    # 2. you also need to include protein-coding genes only

    de_gene_mod = de_gene_df.set_index('hgnc_symbol')
    pgenmi_dict = collections.OrderedDict()
    for de_gene in list(de_gene_df['hgnc_symbol']):
        if not (de_gene in protein_coding_genes):
            continue
        arr = np.zeros(len(tf_dict.keys())+2)
        #print(de_gene_mod.loc[de_gene, 'pvalue'])
        arr[0] = de_gene_mod.loc[de_gene, 'pvalue']
        if math.isnan(de_gene_mod.loc[de_gene, 'pvalue']):
            arr[0] = 1
        arr[1] = 1
        if de_gene in gene_dict.keys():
            for tf in gene_dict[de_gene]:
                ind = tf_dict[tf]
                arr[ind+1] = 1
        pgenmi_dict[de_gene] = arr
        
    # ref for pd.dataframe slicing and indexing: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html
    col_labels = ['PVAL','CONST']+list(tf_dict.keys())
    pgenmi_input = pd.DataFrame.from_dict(pgenmi_dict, orient='index', columns= col_labels)
    pgenmi_input.to_csv(f'{out_path}/tfbs_only_{dist}', sep= '\t', index_label= False)
    
    
#########################################
### TFBS-DiffMark, TFBS-DiffMark-Aggr ###
#########################################

marks = ['K27ac', 'K27me3', 'K4me1', 'K4me3']

# get TF names using (dist, mark) pair having all the TFs (20 TFs)
dist= '1Mb'
mark= 'K27ac'
TFs= []
tfbs_diffmark_dir= f'binary_feature_diffmark/{dist}/{mark}'
for TF in os.listdir(tfbs_diffmark_dir):
    if not '._' in TF:
        TFs.append(TF)
        
df_ref = pd.read_csv(f'{tfbs_diffmark_dir}/{TFs[0]}', sep= '\t')
df_ref.columns = ['PVAL', 'CONST', 'A1', 'A2']
gene_no = df_ref.shape[0]
gene_index = df_ref.index.tolist()


# printing the TF, mark pairs with all-zero evidence (no evidence)
out_path_diffmark = f'{evid_rootdir}/TFBS_DiffMark/'
out_path_diffmark_aggr = f'{evid_rootdir}/TFBS_DiffMarkAggr/'

os.makedirs(out_path_diffmark, exist_ok=True)
os.makedirs(out_path_diffmark_aggr, exist_ok=True)
    
file= open(f'{out_path_diffmark}/missingEvid', 'w')
for dist in dists:
    first_time= False
    i= 0 
    for TF in TFs:
        for mark in marks:
            path = f'binary_feature_diffmark/{dist}/{mark}'
            if TF in os.listdir(path):
                tfbs_diffmark_dir = f'{path}/{TF}'
                df_tf = pd.read_csv(tfbs_diffmark_dir, sep= '\t')
                df_tf.columns= ['PVAL', 'CONST', f'{TF}_{mark}_up', f'{TF}_{mark}_down']
                if not first_time:
                    df_all = df_tf
                    first_time = True
                else:
                    evid_cols = df_tf.columns[2:]
                    df_all = pd.concat([df_all, df_tf[evid_cols]], axis=1)
            else:
                i= i+1
                file.write(f'{mark}, {TF}\n')
                print(mark, TF)
                arr_zeros = np.zeros(shape=(gene_no,2), dtype= int)
                df_zeros = pd.DataFrame(arr_zeros, columns=[f'{TF}_{mark}_up',f'{TF}_{mark}_down'], 
                                        index= gene_index)
                if not first_time:
                    evid_cols = df_ref.columns[:2] # get pvalue and constant
                    df_all = pd.concat([df_ref[evid_cols], df_zeros] ,axis=1)
                    first_time = True
                else:
                    df_all = pd.concat([df_all, df_zeros], axis=1)
                    
    df_all.to_csv(f'{out_path_diffmark}/tfbs_diffMark_{dist}', sep='\t', index_label= False)
    print(dist, i)
    file.write(f'{dist}, {i}\n')
    
    ######################
    # TFBS_DiffMark_aggr #
    ######################
    # generate TFBS overlapping with any mark changing for up/down seperatley evidence
    
    df_aggr= df_all[['PVAL', 'CONST']]
    for TF in TFs:
        tf_up_lbl=[]
        tf_down_lbl=[]
        for lbl in df_all.columns.tolist():
            if TF in lbl and 'up' in lbl:
                if TF == 'EZH2' and (not 'EZH2phosphoT487' in lbl):
                    tf_up_lbl.append(lbl)
                elif TF != 'EZH2':
                    tf_up_lbl.append(lbl)
            if TF in lbl and 'down' in lbl:
                if TF == 'EZH2' and (not 'EZH2phosphoT487' in lbl):
                    tf_down_lbl.append(lbl)
                elif TF != 'EZH2':
                    tf_down_lbl.append(lbl)

        df_tf_up = df_all[tf_up_lbl]
        df_tf_down = df_all[tf_down_lbl]

        # convert series to dataframe
        df_up_aggr = df_tf_up.sum(axis=1).to_frame()
        df_up_aggr.columns = [f'{TF}_up']
        df_down_aggr = df_tf_down.sum(axis=1).to_frame()
        df_down_aggr.columns = [f'{TF}_down']

        df_aggr = pd.concat([df_aggr, df_up_aggr], axis=1)
        df_aggr = pd.concat([df_aggr, df_down_aggr], axis=1)
    evid_cols = df_aggr.columns.tolist()[2:]
    df_aggr[df_aggr[evid_cols] > 0] = 1
    df_aggr.to_csv(f'{out_path_diffmark_aggr}/tfbs_diffMark_aggr_{dist}', sep='\t', index_label= False)
          
file.close()


#########################################
###         TFBS-PresMark             ###
#########################################

marks = ['K27ac', 'K27me3', 'K4me1', 'K4me3']

# get TF names using (dist, mark) pair having all the TFs (20 TFs)
dist= '1Mb'
mark= 'K27ac'
TFs= []
tfbs_diffmark_dir= f'/binary_feature_presmark/{dist}/{mark}'
for TF in os.listdir(tfbs_diffmark_dir):
    if not '._' in TF:
        TFs.append(TF)
        
df_ref= pd.read_csv(f'{tfbs_diffmark_dir}/{TFs[0]}', sep= '\t')
df_ref.columns= ['PVAL', 'CONST', 'A1']
gene_no= df_ref.shape[0]
gene_index= df_ref.index.tolist()

# printing the TF, mark pairs with all-zero evidence (no evidence)
out_path= f'{evid_rootdir}/TFBS_PresMark/'
os.makedirs(out_path, exist_ok=True)
file= open(f'{out_path}/missingEvid', 'w')

for dist in dists:
    first_time= False
    i= 0 
    for TF in TFs:
        for mark in marks:
            path= f'binary_feature_presmark/{dist}/{mark}'
            if TF in os.listdir(path):
                tfbs_diffmark_dir= f'{path}/{TF}'
                df_tf= pd.read_csv(tfbs_diffmark_dir, sep= '\t')
                df_tf.columns= ['PVAL', 'CONST', f'{TF}_{mark}_pres']
                if not first_time:
                    df_all= df_tf
                    first_time= True
                else:
                    evid_cols= df_tf.columns[2:]
                    df_all = pd.concat([df_all, df_tf[evid_cols]], axis=1)
            else:
                i= i+1
                file.write(f'{mark}, {TF}\n')
                print(mark, TF)
                arr_zeros= np.zeros(shape=(gene_no,1), dtype= int)
                df_zeros = pd.DataFrame(arr_zeros, columns=[f'{TF}_{mark}_pres'], 
                                        index= gene_index)
                if not first_time:
                    evid_cols= df_ref.columns[:2] # get pvalue and constant
                    df_all= pd.concat([df_ref[evid_cols], df_zeros] ,axis=1)
                    first_time= True
                else:
                    df_all = pd.concat([df_all, df_zeros], axis=1)
    df_all.to_csv(f'{out_path}/tfbs_presMark_{dist}', sep='\t', index_label= False)
    print(dist, i)
    file.write(f'{dist}, {i}\n')
file.close()


########################################
###           TFBS-DiffAcc           ###
########################################

marks = ['ACC']
# get TF names using (dist, mark) pair having all the TFs (20 TFs)
dist= '1Mb'
mark= 'K27ac'
tfbs_diffmark_dir= f'binary_feature_diffmark/{dist}/{mark}'
TFs= []
for TF in os.listdir(tfbs_diffmark_dir):
    if not '._' in TF:
        TFs.append(TF)
        
tfbs_diffacc_dir= f'binary_feature_diffmrak/{dist}/ACC'
df_ref= pd.read_csv(f'{tfbs_diffacc_dir}/JUND', sep= '\t')
df_ref.columns= ['PVAL', 'CONST', 'A1', 'A2']
gene_no= df_ref.shape[0]
gene_index= df_ref.index.tolist()

# printing the TF, mark pairs with all-zero evidence (no evidence)
out_path= f'{evid_rootdir}/TFBS_DiffACC/'
os.makedirs(out_path, exist_ok=True)
file= open(f'{out_path}/missingEvid', 'w')

for dist in dists:
    first_time= False
    i= 0 
    for TF in TFs:
        path= f'binary_feature_diffmark/{dist}/ACC'
        if TF in os.listdir(path):
            tfbs_diffacc_dir= f'{path}/{TF}'
            df_tf= pd.read_csv(tfbs_diffacc_dir, sep= '\t')
            df_tf.columns= ['PVAL', 'CONST', f'{TF}_ACC_up', f'{TF}_ACC_down']
            if not first_time:
                df_all= df_tf
                first_time= True
            else:
                evid_cols= df_tf.columns[2:]
                df_all = pd.concat([df_all, df_tf[evid_cols]], axis=1)
        else:
            i= i+1
            file.write(f'{TF}\n')
            print(TF)
            arr_zeros= np.zeros(shape=(gene_no,2), dtype= int)
            df_zeros = pd.DataFrame(arr_zeros, columns=[f'{TF}_ACC_up',f'{TF}_ACC_down'], 
                                    index= gene_index)
            if not first_time:
                evid_cols= df_ref.columns[:2] # get pvalue and constant
                df_all= pd.concat([df_ref[evid_cols], df_zeros] ,axis=1)
                first_time= True
            else:
                df_all = pd.concat([df_all, df_zeros], axis=1)
    df_all.to_csv(f'{out_path}/tfbs_diffacc_{dist}', sep='\t', index_label= False)
    print(dist, i)
    file.write(f'{dist}, {i}\n')
file.close()


########################################
###           TFBS-PresAcc           ###
########################################

# generating pgenmi input files for the overlap of TFBS with an accessibility peak in either of the stages
# get TF names using (dist, mark) pair having all the TFs (20 TFs) => comments have not been 
# modified for accessibility data 

dist= '1Mb'
mark= 'ACC'
TFs= []
tfbs_diffmark_dir= f'binary_feature_presmark/{dist}/{mark}'
for TF in os.listdir(tfbs_diffmark_dir):
    if not '._' in TF:
        TFs.append(TF)
        
df_ref= pd.read_csv(f'{tfbs_diffmark_dir}/{TFs[0]}', sep= '\t')
df_ref.columns= ['PVAL', 'CONST', 'A1']
gene_no= df_ref.shape[0]
gene_index= df_ref.index.tolist()

# printing the TF, mark pairs with all-zero evidence (no evidence)
out_path= f'{evid_rootdir}/TFBS_PresACC/'
os.makedirs(out_path, exist_ok=True)
file= open(f'{out_path}/missingEvid', 'w')

for dist in dists:
    first_time= False
    i= 0 
    for TF in TFs:
        path= f'binary_feature_presmark/{dist}/ACC'
        if TF in os.listdir(path):
            tfbs_diffmark_dir= f'{path}/{TF}'
            df_tf= pd.read_csv(tfbs_diffmark_dir, sep= '\t')
            df_tf.columns= ['PVAL', 'CONST', f'{TF}_ACC_pres']
            if not first_time:
                df_all= df_tf
                first_time= True
            else:
                evid_cols= df_tf.columns[2:]
                df_all = pd.concat([df_all, df_tf[evid_cols]], axis=1)
        else:
            i= i+1
            file.write(f'ACC, {TF}\n')
            print(TF)
            arr_zeros= np.zeros(shape=(gene_no,1), dtype= int)
            df_zeros = pd.DataFrame(arr_zeros, columns=[f'{TF}_ACC_pres'], 
                                    index= gene_index)
            if not first_time:
                evid_cols= df_ref.columns[:2] # get pvalue and constant
                df_all= pd.concat([df_ref[evid_cols], df_zeros] ,axis=1)
                first_time= True
            else:
                df_all = pd.concat([df_all, df_zeros], axis=1)                
    df_all.to_csv(f'{out_path}/tfbs_presACC_{dist}', sep='\t', index_label= False)
    print(dist, i)
    file.write(f'{dist}, {i}\n')
file.close()

##########################################################
###    UP-regulated vs DOWN-regulated Gene Separation  ###
###         for up-analysis vs down-analysis           ###
##########################################################

def up_down_filter(df_all, final_de_up_gene, final_de_down_gene):
    df_all_up= df_all.copy(deep= True)
    df_all_down= df_all.copy(deep= True)
    
    df_all_up.PVAL= df_all_up.PVAL/2
    df_all_down.PVAL= df_all_down.PVAL/2
    
    # up-regaulated pgenmi input generation
    # making down-regulated genes insignificant
    df_all_up.loc[final_de_down_gene, 'PVAL']= 1- df_all_up.loc[final_de_down_gene, 'PVAL']
    # down-regaulated pgenmi input generation
    # making up-ragulated genes insignificant
    df_all_down.loc[final_de_up_gene, 'PVAL']= 1- df_all_down.loc[final_de_up_gene, 'PVAL']
    
    return df_all_down, df_all_up

paths= ['TFBS_only/tfbs_only_', 'TFBS_DiffMark/tfbs_diffMark_',
        'TFBS_DiffMarkAggr/tfbs_diffMark_aggr_', 'TFBS_PresMark/tfbs_presMark_',
        'TFBS_DiffACC/tfbs_diffacc_', 'TFBS_PresACC/tfbs_presACC_']

for dist in dists:
    for path in paths:
        path_up= path.split('/')[0]+f'/{dist}/up/'
        path_down= path.split('/')[0]+f'/{dist}/down/'
        
        os.makedirs(path_up, exist_ok=True)
        os.makedirs(path_down, exist_ok=True)        
        print(path_up, path_down)

        df_all= pd.read_csv(f'{path}{dist}', sep= '\t')
        de_gene_down_lst= de_gene_df.loc[de_gene_df['log2FoldChange'] < 0].hgnc_symbol.tolist()
        de_gene_up_lst= de_gene_df.loc[de_gene_df['log2FoldChange'] > 0].hgnc_symbol.tolist()
        final_genes= df_all.index.tolist()
        # find the intersection of up/down with final genes
        final_de_up_gene= sorted(list(set(final_genes) & set(de_gene_up_lst)), reverse= True)
        final_de_down_gene= sorted(list(set(final_genes) & set(de_gene_down_lst)), reverse= True)
        
        df_all_down, df_all_up= up_down_filter(df_all, final_de_up_gene, final_de_down_gene)
        
        df_all_up.to_csv(f'{path_up}/H1_up', sep= '\t', index_label= False)
        df_all_down.to_csv(f'{path_down}/H1_down', sep= '\t', index_label= False)
        
        df_all_up[['PVAL', 'CONST']].to_csv(f'{path_up}/H0_up', sep= '\t', index_label= False)
        df_all_down[['PVAL', 'CONST']].to_csv(f'{path_down}/H0_down', sep= '\t', index_label= False)
