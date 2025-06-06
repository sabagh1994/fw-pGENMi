#!/bin/bash


cd intersection
declare -A dists=( ["1Mb"]=1000000 ["200Kb"]=200000 ["50Kb"]=50000 ["10Kb"]=10000 )

for distname in "${!dists[@]}"
do
	# Get closest gene within e.g. distname 200KB
	distval=${dists[$distname]}
	awk -F $'\t' -v thr=$distval 'BEGIN{OFS=FS} {if ($14 < thr) print $1,$2,$3,$4"|"$13,$5,$6,$7,$8,$9}' coords/tfs_intergenic_protein_coding_closest_gene_coords.tsv >> tfs_intergenic_protein_coding_with_closest_gene_$distname.bed10
	sort -k1,1 -k2,2n tfs_intergenic_protein_coding_with_closest_gene_$distname.bed10 -o tfs_intergenic_protein_coding_with_closest_gene_$distname.bed10.sorted
	# Concatenate intronic and intergenic bed10 files into final file
	cat tfs_intergenic_protein_coding_with_closest_gene_$distname.bed10.sorted > ../tfs_protein_coding_$distname.bed10
	cat tfs_intron_protein_coding.bed10.sorted >> ../tfs_protein_coding_$distname.bed10
	# Format: chr1	856186	857068	SP1|AL645608.1	46	.	140.49	-1	-1
	sort -k1,1 -k2,2n ../tfs_protein_coding_$distname.bed10 -o ../tfs_protein_coding_$distname.bed10.sorted
done
