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
3. **Overlap.** `cd 03_overlap/encode_gencode`, then run `./runme1.sh` followed by `./runme2.sh`. Identifies TF binding sites overlapping with intergenic and intronic regions. Then associates TF binding sites with genes at varying regulatory distances 10Kb, 50Kb, 200Kb and 1Mb.
After running, the content of `03_overlap/encode_gencode` should include the following files, where {dist} represents the regulatory distance.
    ```
    03_overlap/encode_gencode/
    ├── coords
    │   ├── tfs_intergenic_protein_coding_bed_coords.tsv
    │   └── tfs_intron_protein_coding_bed_coords.tsv
    │
    ├── intersection
    │   ├── coords
    │   │   └── tfs_intergenic_protein_coding_closest_gene_coords.tsv
    │   ├── tfs_intergenic_protein_coding.bed10
    │   ├── tfs_intergenic_protein_coding.bed10.sorted
    │   ├── tfs_intergenic_protein_coding_with_closest_gene_{dist}.bed10
    │   ├── tfs_intergenic_protein_coding_with_closest_gene_{dist}.bed10.sorted
    │   ├── tfs_intron_protein_coding.bed10
    │   └── tfs_intron_protein_coding.bed10.sorted
    │ 
    ├── tfs_protein_coding_{dist}.bed10
    ├── tfs_protein_coding_{dist}.bed10.sorted
    
    ```

4. **Epigenetic Marks.** `cd 04_epi_analysis/01_bedtools_mark_acc` then run `./runme.sh` which,\
   **(i).** intersects the replicate measurements for each histone mark/accessibility peaks for each stage of metastasis ($P_0$, noninvasive - $P_6$, metastatic). The directory structure should look like the following after running,

    ```
    04_epi_analysis/01_bedtools_mark_acc/
    ├── final
    │   ├── bed5
    │   │   └── all
    │   │       ├── K27ac.bed (similar for other histone marks and accessibility)
    │   │       
    │   └── narrowPeak
    │       └── all
    │           ├── K27ac.narrowPeak (similar for other histone marks and accessibility)
    |
    ├── K27ac (similar structure for other histone marks and accessibility)
    │   ├── final
    │   │   ├── de_p0_p6.final.narrowPeak
    │   │   ├── de_p0_p6.final.ranked.peak.narrowPeak
    │   │   ├── p0.exclusive.narrowPeak
    │   │   └── p6.exclusive.narrowPeak
    │   ├── p0
    │   │   ├── p0_intersection.final
    │   │
    │   ├── p6
    │   │   ├── p6_intersection.final
    │   │
    │   └── union
    │       └── p0_p6.union.NarrowPeak
    └── Union
        ├── K27ac.bed (similar for other histone marks and accessibility)
    
    ```

   
   **(ii).** Finding the differential mark peaks locations, i.e., gain or loss of peak in transitioning from $P_0$ to $P_6$ (diffmark), as well as the
      presence of mark in either of the stages (presmark). Here "mark" stands for both histone mark and accessibility peaks.
   
5. **Intersect Mark and TF Binding Site.** `cd 04_epi_analysis/02_intersect_tfbs` then run `./runme.sh`. Intersection of **Step4** (diffmark, presmark) with TF binsing sites computed at varying regulatory distances from **Step 3** followed by binarizing the evidence. The directory structure after running looks like the following,
    ```
    04_epi_analysis/02_intersect_tfbs/
    ├── intersections_diffmark
    │   ├── 10Kb (similar for other regulatory distances)
    │       ├── by_mark_10Kb
    │       │   └── binary
    │       │       ├── K27ac
    │       │           ├── ATF3 (similar for other TFs)
    │       │ 
    │       ├── coords_10Kb
    │       │   ├── K27ac.bed (similar for other marks)
    │       │ 
    │       └── gene_sets_10Kb
    │           ├── K27ac.bed (similar for other marks)
    │           └── max_by_tf
    │               ├── binary
    │               │   ├── K27ac.bed (similar for other marks)
    │               │   
    │               └── continuous
    │                   ├── K27ac.bed (similar for other marks)
    │   
    └── intersections_presmark (similar to intersections_diffmark)
        ├── 10Kb
    
    ```

6. **Input Evidence File.** `cd 05_inputgen` then run `runme.sh`. Aggregate the evidence over all TFs, then combines it with DE gene p-values - adjusted for the direction of analysis - to generate the final input files per evidence type $\in$ {DiffMark, DiffMarkAggr, DiffAcc, TFBS-only, PresMark and PresAcc}.
