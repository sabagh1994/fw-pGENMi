## Preprocessing of the Histone Marks and Accessibility Peaks
Follow the instructions below to generate input evidence files by intersecting TF binding sites with histone mark and accesibility peaks.
It is assumed that peak calling has been performed on the raw data.

### Summary
Make sure the following steps have been done before preprocessing.
1. Install [bedtools](https://github.com/sabagh1994/fw-pGENMi/tree/master/README.md#make-venv)
2. Follow the guide in **[Step 3. Download the Data](https://github.com/sabagh1994/fw-pGENMi/blob/master/README.md#download-data)** to download the multi-omics data. The downloaded data will be stored at `preprocess/data` directory. The structure of this directory is represented partially below.
    ```
    preprocess/data/
    ├── ATAcSeq
    │   ├── M0_1101
    │   │   ├── M0_1101.#_hg19_q30_macs_peaks.narrowPeak
    │   │
    │   ├── M6_1101
    │       ├── M6_1101.#_hg19_q30_macs_peaks.narrowPeak
    │
    ├── deseq2
    │   └── p0vsp6_DESeq_processed
    ├── K27ac
    │   ├── Rep1
    │   │   ├── P6-#_macs2_peaks.encodePeak
    │   │   └── SW480-#_macs2_peaks.encodePeak
    │   └── Rep2
    ├── K27me3
    │   ├── Rep1
    ```
3. Run the preprocessing as explained below. 

### Preprocessing

1. `01_gencode`: Gene locations in genome. you should download the right version for human genome.
2. `02_encode`: Download TF Binding profiles from ENCODE and process them.
3. `03_overlap`: Finding the TFBS associated with genes at varying regulatory distances 10Kb, 50Kb, 200Kb and 1Mb.
4. `04_epi_analysis/01_bedtools_mark_acc`: Finding the differential mark peaks locations, i.e., gain/loss of peak in p0/p6 (diffmark), 
                                           as well as the presence of mark in either of the stages noninvasive-p0/metastatic-p6 (presmark).
                                           Here 'mark' stands for both histone mark and accessibility.
                                           
5. `04_epi_analysis/02_intersect_tfbs`: Intersection of step4 (diffmark, presmark) with TFBS computed at varying regulatory distances
                                        from step 3 followed by binarizing the evidence.
6. `05_inputgen`: generating the final input files by aggregating the evidence over all TFs per evidence type, DiffMark, DiffMarkAggr, DiffAcc, TFBS-only, PresMark and PresAcc.
