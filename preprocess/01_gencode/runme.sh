#!/bin/sh

# Goto https://www.gencodegenes.org/#
# From GENCODE website, however over data/human and select GRCh37-mapped releases
# Select the newest version in GENCODE release column
# Copy link from Comprehensive gene annotation for GTF file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz
gunzip gencode.v27lift37.annotation.gtf.gz 

# Following this guide loosely to get intronic and intergenic regions via bedtools
# url: http://crazyhottommy.blogspot.com/2013/05/find-exons-introns-and-intergenic.html

# Extract protein_coding
# Output Format: chr1	HAVANA	gene	65419	71585	.	+	.	gene_id "ENSG00000186092.5_3"; gene_type "protein_coding"; gene_name "OR4F5"; level 2; havana_gene "OTTHUMG00000001094.3_3"; remap_status "full_contig"; remap_num_mappings 1; remap_target_status "overlap";
cat gencode.v27lift37.annotation.gtf | grep "protein_coding" > gencode.v27lift37.protein_coding.gtf

# Get gene bed information for both protein coding and all genes (your choice)
# Output Format: chr1	11869	14409	DDX11L1
./code/parse_gtf.pl gencode.v27lift37.annotation.gtf gene gene_name > gencode.v27lift37.annotation.gene.bed4
./code/parse_gtf.pl gencode.v27lift37.protein_coding.gtf gene gene_name > gencode.v27lift37.protein_coding.gene.bed4

# Get exon bed information for both protein coding and all genes (your choice)
# Output Format: chr1	11869	12227	DDX11L1
./code/parse_gtf.pl gencode.v27lift37.annotation.gtf exon gene_name > gencode.v27lift37.annotation.exon.bed4 
./code/parse_gtf.pl gencode.v27lift37.protein_coding.gtf exon gene_name > gencode.v27lift37.protein_coding.exon.bed4

# Remove GL chr from annotation.gtf for exons and genes
# Output Format: chr1	11869	12227	DDX11L1
cat gencode.v27lift37.annotation.exon.bed4 | egrep -v "^GL" > blah && mv blah gencode.v27lift37.annotation.exon.bed4
cat gencode.v27lift37.annotation.gene.bed4 | egrep -v "^GL" > blah && mv blah gencode.v27lift37.annotation.gene.bed4

# Sort gene bed4 files
for i in *.bed4; do 
	sort -k1,1 -k2,2n $i -o $i.sorted
done

# Merge overlapping exons for sorted bed4 exons files
# Output Format: chr1	11869	12227
mergeBed -i gencode.v27lift37.annotation.exon.bed4.sorted > gencode.v27lift37.annotation.exon.merged.bed4
mergeBed -i gencode.v27lift37.protein_coding.exon.bed4.sorted > gencode.v27lift37.protein_coding.exon.merged.bed4

# Sort gene bed4 files again (you can sort each indep if you want, this re-sorts existing)
for i in *.bed4; do 
	sort -k1,1 -k2,2n $i -o $i.sorted
done

# Define intronic regions as diff of gene and exon.merged bed4 regions
# Output Format: chr1	65433	65520	OR4F5
subtractBed -a gencode.v27lift37.protein_coding.gene.bed4.sorted -b gencode.v27lift37.protein_coding.exon.merged.bed4.sorted > gencode.v27lift37.protein_coding.intron.bed4
subtractBed -a gencode.v27lift37.annotation.gene.bed4.sorted -b gencode.v27lift37.annotation.exon.merged.bed4.sorted > gencode.v27lift37.annotation.intron.bed4 

# Get fetchChromSizes from UCSC toolds
# wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
# put it in code/
# Get lengths for each chromosome in hg19
# Sort chromosome file (strange but it has to happen)
./code/fetchChromSizes hg19 > hg19_chrom_sizes.tsv
sort -k1,1 hg19_chrom_sizes.tsv -o hg19_chrom_sizes.tsv 

# Get intergenic regions
# Output Format: chr1	31109	34554
complementBed -i gencode.v27lift37.annotation.gene.bed4.sorted -g hg19_chrom_sizes.tsv > gencode.v27lift37.annotation.intergenic.bed4
complementBed -i gencode.v27lift37.protein_coding.gene.bed4.sorted -g hg19_chrom_sizes.tsv > gencode.v27lift37.protein_coding.intergenic.bed4

# Sort intergenic bed4 files
for i in *.intergenic.bed4; do
	sort -k1,1 -k2,2n $i -o $i.sorted
done

# Get protein_coding genes
# Output Format: each gene name in a line
cut -f 4 gencode.v27lift37.protein_coding.gene.bed4 | sort -u > protein_coding_genes.list
