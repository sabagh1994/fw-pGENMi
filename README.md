# An Integrated Multi-Omics Approach to Identify Regulatory Mechanisms in Cancer Metastatic Processes

This repository contains the code used for ["An integrated multi-omics approach to identify regulatory mechanisms in cancer metastatic processes"](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02213-x), published in Genome Biology 2021.
## Summary
To get familiar with the study design and be able to reproduce it on a different multi-omics dataset, the following steps should be performed in the given order,
1) Clone the repository to create a local copy of it on your machine.
2) Make a virtual environment and install the required packages.
3) Download the raw data for this study containing the histone mark peaks, ATACseq and differential gene expression statistics from google drive.
   The data can also be accessed from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA659546).
4) Preprocess the raw data, downloaded in step 2, to create input files with different evidence types for pgenmi and fwpgenmi training. 
   Navigate to `preprocess` folder and follow the instructions.
5) Run bash scripts in `exc` folder to perform,
   (i). Cross-validation on pgenmi/fwpgenmi,
   (ii). Train final models using the best hyper parameters setting (regulatory distance, $L_2$ regularization coefficient)
   (iii). Rank the transcription factors for their association to phenotype, e.g., cancer progression in this study.
   (iv). Derive gene sets mediating the impact of highly-ranked TFs on phenotype. These sets can be used as gene signatures in survival analysis.
   (v). Cluster multi-omics profiles of TCGA COAD patients and perform survival analysis using the gene signatures derived from (iv).

Note: The evidence input files used in this study are available on google drive (download by running `exc/01_data_download.sh`), so you can skip preprocessing the raw data (Step 4).


## Step 1. Clone this Repositiry
1. Clone the repository, `git clone https://github.com/sabagh1994/fw-pGENMi.git`
2. Navigate to the cloned directory, `cd fw-pGENMi`

## Step 2. Make a Virtual Environemnt
You need to create your own environment with the required packages installed. There are two ways to create a virtual environment, (1) mamba and (2) python venv module, as described in the following,
1. To make a virtual environment using mamba run `make mamba`. This will download micromamba and install Python 3.11 followed by all the packages 
   required for this project listed in `requirements.txt`.
2. If instead you would like to use python venv, first you need to have python installed. Then run `make venv`.

Note: Bedtools is required for preprocessing. To install it run `make bedtools`.
To activate the virtual environments run `./activate venv` or `./activate mamba` depending on the way the environment was created.

## Step 3. Download the Data


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
