#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


###################
# Command line args
###################

#sqtl viz input file
## File format:
##### 1: variant id
##### 2: cluster id
##### 3: junction start pos
##### 3: junction end pos
##### 4: gene-id
##### 5: Tissue-id
##### 6: Path to genotype file
sqtl_viz_input_file="$1"

# Where to save output
output_dir="$2"


#################
# Reference data
#################
# File containing annotated exons
exon_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/gencode_v26_exons.txt"

# Directory containing LeafCutter splicing counts
# One file for each tissue
leafcutter_dir="/work-zfs/abattle4/lab_data/GTEx_v8_LeafCutter_counts/GTEx_Analysis_v8_sQTL_leafcutter_counts/"


# Loop through sQTLs
while IFS=$'\t' read -r -a myArray
do
	# Extract relevent fields from this sQTL
	snp_id=${myArray[0]}
	cluster_id=${myArray[1]}
	ensamble_id=${myArray[4]}
	tissue_name=${myArray[5]}
	genotype_file=${myArray[6]}

	cluster_to_plot_file=$output_dir"sqtl_viz_"$cluster_id"_"$ensamble_id"_"$snp_id"_"$tissue_name".txt"
	junction_file=$leafcutter_dir$tissue_name"_perind_numers.counts.gz"
	python make_clusters_to_plot_file.py $ensamble_id $cluster_id $genotype_file $junction_file $exon_file $cluster_to_plot_file

	python extract_cluster_counts.py $cluster_to_plot_file $tissue_name $junction_file $output_dir
	
	Rscript visualize_cluster_distribution.R $cluster_to_plot_file $tissue_name $exon_file $output_dir

done <$sqtl_viz_input_file






