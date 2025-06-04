#!/bin/sh
#
# We will be identifying overlapping TFBS with intergenic and intronic regions
# Make soft links for convenience
ln -s ../../02_encode/tfs.bed10.sorted .
ln -s ../../01_gencode/gencode.v27lift37.protein_coding.gene.bed4.sorted .
ln -s ../../01_gencode/gencode.v27lift37.protein_coding.intron.bed4.sorted .
ln -s ../../01_gencode/gencode.v27lift37.protein_coding.intergenic.bed4.sorted .

# Make dir coords (which contain overlap coords for both files) and intersection (which is just the bed file)
mkdir coords intersection
# Output Format: chr1	860124	860604	EZH2phosphoT487	659	.	31.65359	-1.00000	2.28103	chr1	860328	860530	SAMD11
intersectBed -a tfs.bed10.sorted -b gencode.v27lift37.protein_coding.intron.bed4.sorted -wa -wb > coords/tfs_intron_protein_coding_bed_coords.tsv
intersectBed -a tfs.bed10.sorted -b gencode.v27lift37.protein_coding.intergenic.bed4.sorted -wa -wb > coords/tfs_intergenic_protein_coding_bed_coords.tsv 

# Intersect coords files
# Output Format: chr1	91035	91451	CTCF	612	.	19.55880	-1.00000	5.12429
./code/intersect.pl coords/tfs_intergenic_protein_coding_bed_coords.tsv > intersection/tfs_intergenic_protein_coding.bed10
./code/intersect.pl coords/tfs_intron_protein_coding_bed_coords.tsv > intersection/tfs_intron_protein_coding.bed10

# Sort files
cd intersection
for i in *.bed10; do 
	sort -k1,1 -k2,2n $i -o $i.sorted
done


# Get closest gene for integenic
mkdir coords

# bedtools must be installed, for installation check the Makefile at the project base directory
# Output Format: chr1  91035   91451   CTCF    612     .       19.55880        -1.00000        5.12429 chr1    65419   71585   OR4F5   19451
bedtools closest -a tfs_intergenic_protein_coding.bed10.sorted -b ../gencode.v27lift37.protein_coding.gene.bed4.sorted -d -k 20181 > coords/tfs_intergenic_protein_coding_closest_gene_coords.tsv
