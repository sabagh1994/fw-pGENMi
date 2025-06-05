'''
    In this script the input file for coca analysis is generated.
    ---------------------------------------------------------------
    
    Before running this script make sure that 
     1. mediator gene sets have been generated and
     2. Data for TCGA COAD patients has been downloaded and stored in GDC_TCGA_COAD folder.
    
    For each mediator gene set, 
    The intersection of each omic's genes and the mediator gene set is taken.
    For microRNAs, take the intersection of miRNA targest and the mediator gene set.
    The microRNAs are available in both mirTarbase and TCGA data.
'''

import os
import argparse
import pandas as pd
from io_cfg import results_dir, input_dir
from io_cfg import make_abspath, read_yaml

# -------------------------------- #
#        Get the Arguments         #
# -------------------------------- #
parser = argparse.ArgumentParser()
parser.add_argument('--config_path', type=str, required=True) # e.g., 06_cfg_coca_inputgen.yml
args = parser.parse_args()


# -------------------------------- #
#         Get the Configs          #
# -------------------------------- #
cfg_path = args.config_path
config_dict = read_yaml.load(cfg_path)

survival_addr= config_dict.get('survival_addr', 'GDC_TCGA_COAD/Phenotype_Survival/TCGA-COAD.survival.tsv')
gene_expr= config_dict.get('gene_expr', 'GDC_TCGA_COAD/Gene_Expr/RNA-seq_FPKM/TCGA-COAD.htseq_fpkm.tsv')
# miRNA data:
mirna_addr= config_dict.get('mirna_addr', 'GDC_TCGA_COAD/miRNA/TCGA-COAD.mirna.tsv')
gene2ENS= config_dict.get('gene2ENS', 'GDC_TCGA_COAD/Gene_Expr/RNA-seq_FPKM/gencode.v22.annotation.gene.probeMap')
# mir-target database
target_scan= config_dict.get('target_scan', 'GDC_TCGA_COAD/miRNA/mir_gene_target.csv') # mirTarBase was used
# Somatic Mutation:
soma_addr= config_dict.get('soma_addr', 'GDC_TCGA_COAD/Somatic_MuSEVariantAggr/TCGA-COAD.muse_snv.tsv')
gene_counts= config_dict.get('gene_counts', [10, 20, 50, 70, 100, 200])

if not 'input_subdirs' in config_dict.keys():
    # set to default paths to gene sets generated in previous steps
    input_subdirs= [f'{results_dir}/mediator_genes/AllGenes/AllEvidGenes']
    for i in gene_counts:
        input_subdirs.append(f'{results_dir}/mediator_genes/postOddsRatio/Separate/down/postOddsRatio_top{i}')
        input_subdirs.append(f'{results_dir}/mediator_genes/postOddsRatio/Separate/up/postOddsRatio_top{i}')
        input_subdirs.append(f'{results_dir}/mediator_genes/RatioOfPostOddsRatio/Separate/down/topGene{i}_topTF5')
        input_subdirs.append(f'{results_dir}/mediator_genes/RatioOfPostOddsRatio/Separate/up/topGene{i}_topTF5')

        input_subdirs.append(f'{results_dir}/mediator_genes/postOddsRatio/Union/postOddsRatio_top{i}')
        input_subdirs.append(f'{results_dir}/mediator_genes/RatioOfPostOddsRatio/Union/topGene{i}_topTF5')
        #input_subdirs.append(f'DEGenes/DEGenes_top{i}')
else:
    input_subdirs = config_dict['input_subdirs']
    for i in range(len(input_subdirs)):
        input_subdir = input_subdirs[i]
        input_subdir = make_abspath(input_subdir, results_dir)
        input_subdirs[i] = input_subdir
       
outdir = config_dict.get('outdir', 'coca_results/coca_input')
#input_rootdir= config_dict.get('input_rootdir', './GeneSets')

# ------------------------------------------------------------------------------ #
#                check if relative or abs path is provided                       #
# ------------------------------------------------------------------------------ #
survival_addr = make_abspath(survival_addr, input_dir) 
gene_expr = make_abspath(gene_expr, input_dir)
mirna_addr = make_abspath(mirna_addr, input_dir)
gene2ENS = make_abspath(gene2ENS, input_dir)
target_scan = make_abspath(target_scan, input_dir)
soma_addr = make_abspath(soma_addr, input_dir)
gene_counts = make_abspath(gene_counts, input_dir)

outdir = make_abspath(outdir, results_dir)
# ------------------------------------------------------------------------------ #
    
for input_subdir in input_subdirs:
    # load the gene set:
    #med_gene_addr= f'{input_rootdir}/{input_subdir}'
    med_gene_addr = input_subdir
    df_med_gene= pd.read_csv(med_gene_addr, sep= '\t', header= None)
    genes_original= df_med_gene[0].tolist()

    # create the output directory
    reldir = input_subdir.split(results_dir)[-1]
    #out_dir= f'{outdir}/{input_subdir}/'
    out_dir = f'{outdir}/{reldir}'
    os.makedirs(out_dir, exist_ok= True)
    
    
    df_surv= pd.read_csv(survival_addr, sep= '\t')
    df_expr= pd.read_csv(gene_expr, sep= '\t') # ens id is in one column and patients are column labels
    df_mirna= pd.read_csv(mirna_addr, sep= '\t') # mir id is in one column and patients are column labels

    df_soma= pd.read_csv(soma_addr, sep= '\t')
    df_soma_bin= pd.crosstab(df_soma.gene, df_soma.Sample_ID)
    # make genes as a column not index
    df_soma_bin.reset_index(level=0, inplace=True) # gene names are in one column and patients are column labels

    df_ens= pd.read_csv(gene2ENS, sep= '\t') # ens id and gene names are in two seperate columns
    # Drop all the rows that have any nan value
    df_surv.dropna(inplace= True)
    df_expr.dropna(inplace= True)
    df_mirna.dropna(inplace= True)
    
    # take the intersection of mir in database and the mir in our data
    # then for the mir in intersection look at the target genes

    df_mir_targets= pd.read_csv(target_scan, sep= ',', header= None)
    #df_mir_targets # MIRT000002	hsa-miR-20a-5p	HIF1A
    df_mir_targets[1]= df_mir_targets[1].str.lower()
    sample_mirs= [s.lower() for s in df_mirna['miRNA_ID'].tolist()]
    mirtar_mirs= [s.lower() for s in df_mir_targets[1].tolist()]
    # taking the intersection of mir in the mir-target database and the TCGA mir data available
    com_mir= list(set(sample_mirs) & set(mirtar_mirs))
    
    df_mir_target_final= df_mir_targets[df_mir_targets[1].isin(com_mir)]
    mir_targets= df_mir_target_final[2].tolist() # there may be repeats of the same gene
    
    # intersect the mir target genes with mediator genes 
    # get the mirnas that target the genes in the intersection
    # ----------------------------------------------------------
    mir_genes_final= list(set(mir_targets) & set(genes_original))
    # final mirnas to be used (there might be repeats)
    final_mirnas= df_mir_target_final[df_mir_target_final[2].isin(mir_genes_final)][1].tolist()
    
    # keeping the ensembles associated with the mediator gene set 
    # Attn: some genes are associated with several ensemble ids 
    # (ids.subid) you can only keep one of them or all of them.
    # to keep the first one use df_ens.drop_duplicates(subset=['gene'], keep='first')

    # mapping mdeiator gene names to ensemble ids provided
    final_ens= df_ens[df_ens['gene'].isin(genes_original)]['id'].tolist() # will be used for gene expression and somatic data
    df_expr_final= df_expr[df_expr['Ensembl_ID'].isin(final_ens)]
    df_surv_final= df_surv[['sample', '_OS', '_EVENT']]
    
    patients_expr= df_expr_final.columns.tolist()[1:]
    patients_surv= df_surv_final['sample'].tolist()
    patients_mirna= df_mirna.columns.tolist()[1:]
    patients_soma= df_soma_bin.columns.tolist()[1:]
    
    # intersection of patients
    final_patients= list(set(patients_expr) & set(patients_surv))
    final_patients= list(set(final_patients) & set(patients_mirna))
    final_patients= list(set(final_patients) & set(patients_soma))

    df_expr_final= df_expr_final[['Ensembl_ID'] + final_patients]
    df_surv_final= df_surv_final[df_surv_final['sample'].isin(final_patients)]
    df_mirna_final= df_mirna[['miRNA_ID'] + final_patients]
    # restricting the mirs to the ones that target the mediator genes
    df_mirna_final= df_mirna_final[df_mirna_final['miRNA_ID'].isin(final_mirnas)]
    
    # ref: https://stackoverflow.com/questions/37125174/accessing-every-1st-element-of-pandas-dataframe-column-containing-lists
    # removing the .No from the end of ens name
    df_expr_final['Ensembl_ID']= df_expr_final['Ensembl_ID'].str.split('.').map(lambda x: x[0])

    # saving the gene expression input
    df_expr_final.to_csv(out_dir+'/geneExprXPatient', sep= '\t', index= False)
    df_surv_final.to_csv(out_dir+'/PatientXSurvival', sep= '\t', index= False)
    df_mirna_final.to_csv(out_dir+'/miRExprXPatient', sep= '\t', index= False)
    
    # somatic mutation data:
    # only consider ensembles associated with mediator genes
    df_ens_gene= df_ens[df_ens['gene'].isin(genes_original)] # only the genes that are mediators
    df_soma_bin= df_soma_bin[df_soma_bin.gene.isin(df_ens_gene.gene.tolist())]
    _genes= df_soma_bin.gene
    df_ens_gene= df_ens_gene[df_ens_gene.gene.isin(_genes)]
    # each gene is associated with several ensemble ids
    df_ens_gene= df_ens_gene.drop_duplicates(subset=['gene'], keep= 'first')
    
    # df_soma_bin has gene names not ensemble ids
    df_ens_gene= df_ens_gene.sort_values(by= 'gene')
    df_soma_bin= df_soma_bin.sort_values(by= 'gene')
    # replacing the gene names with ensemble ids in somatic data
    df_soma_bin['gene']= df_ens_gene.id.values
    
    df_soma_bin['gene']= df_soma_bin['gene'].str.split('.').map(lambda x: x[0])
    df_soma_bin= df_soma_bin[['gene']+final_patients]
    df_soma_bin.to_csv(out_dir+'/geneXpatient_somatic', sep= '\t', index= False)
    
    file = open(out_dir+'/log.txt','w')
    file.write(f'#patients\t{len(final_patients)}\n')
    file.write(f'#mediator genes\t{len(genes_original)}\n')
    file.write(f'somatic mutation shape:\t{df_soma_bin.shape[0]},{df_soma_bin.shape[1]}\n')
    file.write(f'miRNA shape:\t{df_mirna_final.shape[0]},{df_mirna_final.shape[1]}\n')
    file.write(f'gene expression shape:\t{df_expr_final.shape[0]},{df_expr_final.shape[1]}\n')
    file.close()
    
