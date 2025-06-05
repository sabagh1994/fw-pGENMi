# An Integrated Multi-Omics Approach to Identify Regulatory Mechanisms in Cancer Metastatic Processes

This repository contains the code used for ["An integrated multi-omics approach to identify regulatory mechanisms in cancer metastatic processes"](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02213-x), published in Genome Biology 2021.
## Summary
To get familiar with the study design and be able to reproduce it on a different multi-omics dataset, the following steps should be performed in the given order,
1) Clone the repository to create a local copy of it on your machine.
2) Make a virtual environment and install the required packages.
3) Download the multi-omics data for this study containing the histone mark peaks, accessibility peaks and differential gene expression statistics from google drive.
   The *raw* data can be accessed from [NCBI Sequence Read Archive (SRA)-PRJNA659546](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA659546).
4) Preprocess the multi-omics data, downloaded from google drive in step 2, to create input files with different evidence types for pgenmi and fwpgenmi training. 
   Navigate to `preprocess` folder and follow the instructions.
5) Run bash scripts in `exc` folder to perform,\
   i. Cross-validation on pgenmi/fwpgenmi.\
   ii. Train final models using the best hyper parameters setting (regulatory distance, $L_2$ regularization coefficient).\
   iii. Rank the transcription factors for their association to phenotype, e.g., cancer progression in this study.\
   iv. Derive gene sets mediating the impact of highly-ranked TFs on phenotype. These sets can be used as gene signatures in survival analysis.\
   v. Cluster multi-omics profiles of TCGA COAD patients and perform survival analysis using the gene signatures derived from (iv).

**Note:** The evidence input files used in this study are available on google drive (download by running `exc/01_data_download.sh`), so you can skip preprocessing the multi-omics data (Step 4).


## Step 1. Clone this Repositiry
1. Clone the repository, `git clone https://github.com/sabagh1994/fw-pGENMi.git`
2. Navigate to the cloned directory, `cd fw-pGENMi`

## Step 2. Make a Virtual Environemnt
You need to create your own environment with the required packages installed. There are two ways to create a virtual environment, (1) mamba and (2) python venv module, as described in the following,
1. To make a virtual environment using mamba run `make mamba`. This will download micromamba and install Python 3.11 followed by all the packages 
   required for this project listed in `requirements.txt`.
2. If instead you would like to use python venv, first you need to have python installed. Then run `make venv`.
3. Bedtools is required for preprocessing. To install it run `make bedtools`.

To activate the virtual environments run `./activate venv` or `./activate mamba` depending on the way the environment was created.

## Step 3. Download the Data
Run `./exc/01_data_download.sh` to download the following data from google drive,
1. Input evidence files to train pgenmi/fwpgenmi models. The evidence files are separated by evidence type, e.g., TFBS_DiffMark,
   regulatory distance, and direction of analysis. The files can be accessed at `./input` directory after download.
2. Multi-omics data used to generate the inputs downloaded in (1). The data includes histone mark profiles, accessibility peaks and differential gene expression analysis for noninvasive (p0) and metastatic stages (p6). The data will be downloaded to `./preprocess/data` directory. Users interested in creating the input evidence files from scratch should download this data.

**Note** The *raw* data has been deposited to [NCBI Sequence Read Archive (SRA)-PRJNA659546](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA659546). But you only need the multi-omics data downloaded from google drive to generate the input evidence files. 

## Step 4. Preprocess the Multi-Omics Data
Navigate to `preprocess` directory by running `cd preprocess`. Then follow the instructions to, 
1. Download TF ChIP-seq profiles from Encode,
2. Intersect the TF bindings sites with epigenetic mark/accessibility peaks, and
3. Create the input evidence files for different evidence types, direction of analysis, and regulatory distances.
Note that performing preprocessing from scratch is not required for the next steps.

## Step 5. Training pgenmi/fwpgenmi Models
Before training the models using the configs in `configs` directory, make sure that the input evidence files are located at `./input` and divided by evidence type. Downloading the data in **Step 3** directly saves the input files to `./input` directory. All the scripts for this step are located at `./exc` directory. These bash scripts run python codes that are configured with config files located at `./configs` directory.\


**(A) Cross validation and TF Rankings.**

Run `./exc/02_cross_valid.sh pgenmi`, or `./exc/02_cross_valid.sh fwpgenmi`. This script performs,
1. cross validation for different values of two hyperparameters, regulatory distance and $L_2$-regularization coefficient. This is done by running `src/cv_pgenmi.py`
   on `configs/01_cfg_cv_{model_type}_H1.yml` and `configs/01_cfg_cv_{model_type}_H0.yml` config files, where model_type $\in$ {pgenmi, fwpgenmi}.
2. Aggregate the cross validation results to pick the best hyperparameter setting, training the model on the best setting, and ranking the TFs.
   This is done by running `scr/cvaggr_tr_evidrank.py` using the config `configs/03_cfg_cvaggr_tr_evidrank_{model_type}.yml`.

**(B) Get the Mediator Genes**

Run `./exc/03_mediator_geneset.sh` to get sets of genes mediating the impact of TFs on the phenotype, e.g., cancer progression. 
Mediator genes are chosen based on post odds ratio (POR) or Ratio of Post Odds Ratio (RPOR). Read the paper for a detailed explanation of these two criteria.
This script runs `src/geneset_gen.py` with the config `configs/05_cfg_mediatorgene.yml` to generate mediator gene sets.

**(C) Clustering and Survival Analysis**

1. First you need to download the multi-omics profiles for TCGA COAD patients from [UCSAC Xena](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Colon%20Cancer%20(COAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443). Once downloaded move all data to the same folder named `GDC_TCGA_COAD`.
The structure of this folder should look like the following. You can either move `GDC_TCGA_COAD` to `./input` directory or provide its relative path to the `./input` directory in the config file used for generating the input for clustering (`configs/06_cfg_coca_inputgen.yml`) as explained below.
```
GDC_TCGA_COAD/
├── CopyNo_segment
│   └── TCGA-COAD.masked_cnv.tsv
├── Gene_Expr
│   ├── RNA-seq_Counts
│   │   ├── gencode.v22.annotation.gene.probeMap
│   │   └── TCGA-COAD.htseq_counts.tsv
│   └── RNA-seq_FPKM
│       ├── gencode.v22.annotation.gene.probeMap
│       └── TCGA-COAD.htseq_fpkm.tsv
├── Methyl_450
│   ├── illuminaMethyl450_hg38_GDC
│   └── TCGA-COAD.methylation450.tsv
├── miRNA
│   ├── hsa_MTI.xlsx
│   ├── mir_gene_target.csv
│   └── TCGA-COAD.mirna.tsv
├── Phenotype_Survival
│   ├── TCGA-COAD.GDC_phenotype.tsv
│   └── TCGA-COAD.survival.tsv
└── Somatic_MuSEVariantAggr
    └── TCGA-COAD.muse_snv.tsv
```

2. Run `src/04_coca.sh` to generate inputs for clustering using different sets of mediator genes, cluster the multi-omics profiles of TCGA COAD patients and perform survival analysis for each cluster of patients. This is done by running the following scripts,\
   **a.** `scr/coca_inputgen.py` using config `configs/06_cfg_coca_inputgen.yml`. Generates different input files by subsetting the multi-omics profiles to different gene sets. For details please read the paper.\
   **b.** `scr/KnowEng.py` using config `results/coca_results/network`. Downloads the biological knowledge network for homo sapiens, which is used for network-guided clustering of somatic mutation profiles. The networks will be stored at `results/coca_results/network`.\
   **c.** `scr/coca_analysis.py` using config `configs/07_cfg_coca_analysis.yml`. Performs clustering of multi-omics profiles subset to varying gene signatures followed by survival analysis on the derived clusters. Different parameters for clustering such as number of clusters and degree of network smoothing can be explored.\
   **d.** `src/coca_stats_aggr.py` using config `configs/08_cfg_coca_aggr.yml`. Aggregate the results of clustering and survival analysis in (c) performed with different settings.

### Input Evidence File Format

### Paths in the Config Files
To run the exact pipeline in **Step 5** on your own input evidence file, the following requirements should be met,
1. The input evidence files should have the same directory architecture as the `input` folder. The directory structure is
   `{evid_rootdir}/{evid_type}/{dist}/{dirc}/H*_{dirc}`, where\
   **a.** evid_rootdir is the root directory to all evidence types.\
   **b.** evid_type is the directpry named by the evidence type $\in$ {TFBS_DiffMark, TFBS_only, TFBS_DiffACC, etc}\
   **c.** dist is the directory named by the regulatory distance $\in$ {10Kb, 50Kb, 200Kb, 1Mb}\
   **d.** dirc is the directory named by the direction of analysis $\in$ {up, down}
   Note that you should always provide absolute path to the `evid_rootdir` or its relative path to the `input` folder in the cloned repo.
   
2. By default the results of all runs will be stored at the `./results` directory. To set a different path, you should update the configs
   with the abosolute path to your desired directory or its relative path to `./results` in the cloned repo.

See the config files for full instructions.
   


<!--- # Data:

To get the processed evidence files used as the input of pgenmi/fwpgenmi
download, wget, ... . This will result in a folder named "input" in the fw-pgenmi folder.

notes on the hierarchy of evidence files input/(evidence_type_name)/(evidence_distance)/(direction_of_analysis)/H0_(direction_of_analysis) or H1_(direction_of_analysis) 

also some instructions on the format of the input (inluding the column tags TF_mark_direction) that the mark order should not change 

# Step 1: Cross Validation:

Cross validation is required to find the best evidence distance and regularization coefficient for each evidence type.
To do cross validation for pgenmi navigate to pgenmi folder (cd pgenmi).
To do cross validation for fwpgenmi navigate to fwpgenmi folder (cd fwpgenmi).

After cross validation the model is trained on the entire dataset which is then used to rank the transcription factors underlying the cancer progression.

# Step 2: Getting mediator genes:

After training the model in step 1, navigate to mediator_genest folder to get the genes mediating the impact of transcription factors on the phenotype, e.g. cancer progression. This can be done based on Posterior Odds Ratio (POR) or Ratio of Posterior Odds Ratios (RPORs).

# Step 3: Survival Analysis:

To cluster multi-omics profiles of TCGA-COAD (colorectal cancer) patients, navigate to coca folder. Most of the functions used for this step are adapted from https://knoweng.github.io.
-->
