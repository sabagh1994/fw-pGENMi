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
Preprocessing should be done in the following order,
1. **Gencode.** `cd 01_gencode` then run `./runme.sh`. Downloads the human genome to get the gene coordinates and produces a list of protein coding genes. Make sure that you download the right version of human geneome for your application.
2. **Encode.** `cd 02_encode` then run `./runme.sh`. Downloads TF Binding profiles from ENCODE for HCT116 cell line. The final processed bed file is `02_encode/tfs.bed10.sorted` with the following format for each line,
   ```
   chr1    91035   91451   CTCF    612     .       19.55880        -1.00000        5.12429 208
   ```
4. **Overlap.** `cd 03_overlap`, then run `./runme1.sh` followed by `./runme2.sh`. Identifies TF binding sites overlapping with intergenic and intronic regions. Then associates TF binding sites with genes at varying regulatory distances 10Kb, 50Kb, 200Kb and 1Mb.
5. **Epigenetic Marks.** `cd 04_epi_analysis/01_bedtools_mark_acc` then run `./runme.sh` which,\
   **i.** intersects the replicate measurements for each histone mark/accessibility peaks for each stage of metastasis (p0, noninvasive - p6, metastatic)\
   **ii.** Finding the differential mark peaks locations, i.e., gain or loss of peak in transitioning from p0 to p6 (diffmark), as well as the
      presence of mark in either of the stages (presmark). Here "mark" stands for both histone mark and accessibility peaks.
   
6. **Intersect Mark and TF Binding Site.** `cd 04_epi_analysis/02_intersect_tfbs` then run `./runme.sh`. Intersection of **Step4** (diffmark, presmark) with TF binsing sites computed at varying regulatory distances from **Step 3** followed by binarizing the evidence.
   
7. `05_inputgen`: generating the final input files by aggregating the evidence over all TFs per evidence type, DiffMark, DiffMarkAggr, DiffAcc, TFBS-only, PresMark and PresAcc.
