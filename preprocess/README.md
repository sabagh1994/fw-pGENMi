### Preprocessing of the Histone Marks and Accessibility Peaks

To create the input for running pgenmi and fwpgenmi, execute the codes in each directory in the following order.
You need to have [bedtools](https://bedtools.readthedocs.io/en/latest/) installed to preprocess the peaks and their overlaps with TF binding sites.
See the instructions in Makefile in the project base directory for more details. 
add bedtools to path for easier access, 'export PATH=$PATH:${PROJBASE}/software/bedtools2/bin # bedtools'

1. `01_gencode`: Gene locations in genome. you should download the right version for human genome.
2. `02_encode`: Download TF Binding profiles from ENCODE and process them.
3. `03_overlap`: Finding the TFBS associated with genes at varying regulatory distances 10Kb, 50Kb, 200Kb and 1Mb.
4. `04_epi_analysis/01_bedtools_mark_acc`: Finding the differential mark peaks locations, i.e., gain/loss of peak in p0/p6 (diffmark), 
                                           as well as the presence of mark in either of the stages noninvasive-p0/metastatic-p6 (presmark).
                                           Here 'mark' stands for both histone mark and accessibility.
                                           
5. `04_epi_analysis/02_intersect_tfbs`: Intersection of step4 (diffmark, presmark) with TFBS computed at varying regulatory distances
                                        from step 3 followed by binarizing the evidence.
6. `05_inputgen`: generating the final input files by aggregating the evidence over all TFs per evidence type, DiffMark, DiffMarkAggr, DiffAcc, TFBS-only, PresMark and PresAcc.

