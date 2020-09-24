import numpy as np 
import os
import sys
import pdb
import gzip




def extract_raw_cluster_jxn_data_structure(tissue_specific_jxn_file):
	# Used to skip header
	head_count = 0
	# Initialize cluster_jxn_data_structure
	cluster_jxn_data_structure = {}
	jxn_names_data_structure = {}
	# Stream input file
	f = gzip.open(tissue_specific_jxn_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			#samples = np.asarray(data[1:])
			samples = []
			for ele in data:
				indi = ele.split('-')[0] + '-' + ele.split('-')[1]
				samples.append(indi)
			samples = np.asarray(samples)
			continue
		# Standard line
		# Extract Jxn info from current line
		key = data[0]
		key_info = key.split(':')
		chrom_num = key_info[0]
		start = key_info[1]  # 5' ss
		end = key_info[2]  # 3' ss
		cluster_id = key_info[3]  # Name of cluster

		jxn_read_counts = np.asarray(data[1:]).astype(float)  # Vector of raw read counts for this junction

		# Add jxn to  cluster_jxn_data_structure
		if cluster_id not in cluster_jxn_data_structure:  # If we've never seen this cluster before
			cluster_jxn_data_structure[cluster_id] = []
			jxn_names_data_structure[cluster_id] = []
		cluster_jxn_data_structure[cluster_id].append(jxn_read_counts)  # Add read counts from current junction
		jxn_name = data[0].split(':')[0] + ':' + data[0].split(':')[1] + ':' + data[0].split(':')[2] + ':' + data[0].split(':')[3]
		jxn_names_data_structure[cluster_id].append(jxn_name)
	f.close()
	# Convert from list of arrays to matrix (in each cluster)
	# Loop through all clusters
	all_clusters = cluster_jxn_data_structure.keys()
	for cluster_id in all_clusters:
		# Create matrix from list of arrays
		jxn_matrix = np.transpose(np.asmatrix(cluster_jxn_data_structure[cluster_id]))
		# Add this new matrix to the data structure
		cluster_jxn_data_structure[cluster_id] = {}
		cluster_jxn_data_structure[cluster_id]['samples'] = samples
		cluster_jxn_data_structure[cluster_id]['jxn_matrix'] = jxn_matrix
		cluster_jxn_data_structure[cluster_id]['jxn_names'] = jxn_names_data_structure[cluster_id]
	return cluster_jxn_data_structure, samples

# Remove clusters with more than $max_number_of_junctions_per_cluster
def max_number_of_jxns_filter_ignore_genes(cluster_jxn_data_structure, max_number_of_junctions_per_cluster):
	# Initialize new cluster_jxn_data_structure
	new_jxn_structure = {}
	# Loop through clusters
	for cluster_id in cluster_jxn_data_structure.keys():
		samples = cluster_jxn_data_structure[cluster_id]['samples']
		jxn_mat = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
		jxn_names = cluster_jxn_data_structure[cluster_id]['jxn_names']
		# Get dimensionality of this cluster
		N, K = jxn_mat.shape
		# If there are more than max_number_of_junctions_per_cluster
		if K > max_number_of_junctions_per_cluster:
			pass
		else: # Less than or equal
			new_jxn_structure[cluster_id] = {}
			new_jxn_structure[cluster_id]['samples'] = samples
			new_jxn_structure[cluster_id]['jxn_matrix'] = jxn_mat
			new_jxn_structure[cluster_id]['jxn_names'] = jxn_names
	return new_jxn_structure

# Remove samples with less than $min_reads_per_sample_in_cluster reads (summed acrosss all junctions in cluster)
def min_reads_per_sample_filter(cluster_jxn_data_structure, min_reads_per_sample_in_cluster):
	# Initialize new data structure (post filtering)
	new_jxn_structure = {}
	# Loop through all clusters
	all_clusters = cluster_jxn_data_structure.keys()
	for cluster_id in all_clusters:
		# Extract relevent fields from cluster
		samples = cluster_jxn_data_structure[cluster_id]['samples']
		jxn_mat = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
		# get dimensionality of cluster
		N, K = jxn_mat.shape
		# Itialize vectors of samples and jxn mat (post filtering)
		new_samples = []
		new_jxn_mat = []
		# Loop through samples
		for n in range(N):
			# Vector spanning this samples junction counts for this cluster
			sample_counts = np.squeeze(np.asarray(jxn_mat[n,:]))
			sample_id = samples[n]
			# Sum of reads across all junctions for this cluster
			total_sample_read_counts = np.sum(sample_counts)
			# Check if sample passes filter
			if total_sample_read_counts >= min_reads_per_sample_in_cluster:
				new_samples.append(sample_id)
				new_jxn_mat.append(sample_counts)
		# Convert from list of arrays to matrices
		new_jxn_mat = np.asmatrix(new_jxn_mat)
		new_samples = np.asarray(new_samples)
		# Add filtered matrices to jxn_structure
		new_jxn_structure[cluster_id] = {}
		new_jxn_structure[cluster_id]['samples'] = new_samples
		new_jxn_structure[cluster_id]['jxn_matrix'] = new_jxn_mat
	return new_jxn_structure


def add_intercept_covariate_to_cluster_jxn_data_structure(cluster_jxn_data_structure):
	# Loop through all clusters
	all_clusters = cluster_jxn_data_structure.keys()
	for cluster_id in all_clusters:
		num_samples = len(cluster_jxn_data_structure[cluster_id]['samples'])
		cov_mat = np.ones((num_samples, 1))
		cluster_jxn_data_structure[cluster_id]['covariate_matrix'] = cov_mat
	return cluster_jxn_data_structure


# Add Covariate Matrix to cluster_jxn_data_structure
def add_covariates_to_cluster_jxn_data_structure(cluster_jxn_data_structure, covariate_method):
	if covariate_method == 'none':
		cluster_jxn_data_structure = add_intercept_covariate_to_cluster_jxn_data_structure(cluster_jxn_data_structure)
	return cluster_jxn_data_structure


def create_cluster_based_data_structure(tissue_specific_jxn_file, max_number_of_junctions_per_cluster):
	#Get raw data structure
	#Also get samples, this is the maximum possible samples after filtering
	cluster_jxn_data_structure, samples = extract_raw_cluster_jxn_data_structure(tissue_specific_jxn_file)

	# Remove clusters with more than $max_number_of_junctions_per_cluster
	cluster_jxn_data_structure = max_number_of_jxns_filter_ignore_genes(cluster_jxn_data_structure, max_number_of_junctions_per_cluster)

	# Remove samples with less than $min_reads_per_sample_in_cluster reads (summed acrosss all junctions in cluster)
	# cluster_jxn_data_structure = min_reads_per_sample_filter(cluster_jxn_data_structure, min_reads_per_sample_in_cluster)

	# Add Covariate Matrix to cluster_jxn_data_structure
	
	return cluster_jxn_data_structure, samples

def extract_cluster_info(clusters_to_plot_file):
	clusters = {}
	head_count = 0
	f = open(clusters_to_plot_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		var_pos = data[3]
		geno_0_indi = data[4].split(',')
		geno_1_indi = data[5].split(',')
		geno_2_indi = data[6].split(',')
		strand = data[7]
		clusters[cluster_id] = {}
		clusters[cluster_id]['geno_0_indi'] = geno_0_indi
		clusters[cluster_id]['geno_1_indi'] = geno_1_indi
		clusters[cluster_id]['geno_2_indi'] = geno_2_indi
		clusters[cluster_id]['variant_position'] = var_pos
		clusters[cluster_id]['strand'] = strand
	f.close()
	return clusters

def save_outlier_jxns_to_output_file(individuals, jxn_names, jxn_matrix, ordered_individuals, output_file):
	if jxn_matrix.shape[0] != len(ordered_individuals):
		print('assumption error!!')
		pdb.set_trace()
	# Filter matrix to include only individuals in individuals
	valid_rows = []
	indis = []
	for i,val in enumerate(ordered_individuals):
		if val in individuals:
			if np.sum(jxn_matrix[i,:]) > 0:
				valid_rows.append(i)
				indis.append(val)
	if len(valid_rows) == 0:
		print('assumption error!')
		pdb.set_trace()
	filtered_matrix = jxn_matrix[valid_rows, :]
	if filtered_matrix.shape[0] != len(indis):
		print('assumption errorro!')
		pdb.set_trace()
	# Start writing to output file
	t = open(output_file,'w')
	# Print header
	t.write('individual\t' + '\t'.join(jxn_names) + '\n')
	for i,val in enumerate(indis):
		t.write(val + '\t' + '\t'.join(np.squeeze(np.asarray(filtered_matrix[i,0:])).astype(int).astype(str)) + '\n')
	for i,val in enumerate(indis):
		t.write(val + '\t' + '\t'.join(np.squeeze(np.asarray(filtered_matrix[i,0:])).astype(int).astype(str)) + '\n')
	t.close()

def save_jxns_to_output_file(geno_0_indis, geno_1_indis, geno_2_indis, jxn_names, jxn_matrix, ordered_individuals, output_file):
	if jxn_matrix.shape[0] != len(ordered_individuals):
		print('assumption error!!')
		pdb.set_trace()
	# Filter matrix to include only individuals in individuals
	geno_0_valid_rows = []
	geno_0_ordered_indis = []
	geno_1_valid_rows = []
	geno_1_ordered_indis = []
	geno_2_valid_rows = []
	geno_2_ordered_indis = []
	for i,val in enumerate(ordered_individuals):
		if val in geno_0_indis:
			if np.sum(jxn_matrix[i,:]) >= 0:
				geno_0_valid_rows.append(i)
				geno_0_ordered_indis.append(val)
		if val in geno_1_indis:
			if np.sum(jxn_matrix[i,:]) >= 0:
				geno_1_valid_rows.append(i)
				geno_1_ordered_indis.append(val)
		if val in geno_2_indis:
			if np.sum(jxn_matrix[i,:]) >= 0:
				geno_2_valid_rows.append(i)
				geno_2_ordered_indis.append(val)
	if len(geno_0_valid_rows) == 0 or len(geno_1_valid_rows) == 0 or len(geno_2_valid_rows) == 0:
		print('assumption error!')
		pdb.set_trace()
	geno_0_filtered_matrix = jxn_matrix[geno_0_valid_rows, :]
	geno_1_filtered_matrix = jxn_matrix[geno_1_valid_rows, :]
	geno_2_filtered_matrix = jxn_matrix[geno_2_valid_rows, :]

	# Start writing to output file
	t = open(output_file,'w')
	# Print header
	t.write('individual\t' + '\t'.join(jxn_names) + '\n')
	for i,val in enumerate(geno_0_ordered_indis):
		t.write(val + '\t' + '\t'.join(np.squeeze(np.asarray(geno_0_filtered_matrix[i,0:])).astype(int).astype(str)) + '\n')
	for i,val in enumerate(geno_1_ordered_indis):
		t.write(val + '\t' + '\t'.join(np.squeeze(np.asarray(geno_1_filtered_matrix[i,0:])).astype(int).astype(str)) + '\n')
	for i,val in enumerate(geno_2_ordered_indis):
		t.write(val + '\t' + '\t'.join(np.squeeze(np.asarray(geno_2_filtered_matrix[i,0:])).astype(int).astype(str)) + '\n')
	t.close()

def revise_jxn_names(jxn_names):
	new = []
	for jxn_name in jxn_names:
		new_name = jxn_name.split(':')[1] + ':' + jxn_name.split(':')[2]
		new.append(new_name)
	return new

clusters_to_plot_file = sys.argv[1]
tissue_name = sys.argv[2]
tissue_specific_junction_file = sys.argv[3]
visualize_cluster_distribution_dir = sys.argv[4]

# NOTE: THERE ARE SOME CLUSTER REPEATS HERE. And we only plot one indi per cluster
cluster_struct = extract_cluster_info(clusters_to_plot_file)

# Extract junction data and place in compact data structure
# Keys are cluster_ids and values are jxn_counts
cluster_jxn_data_structure, all_samples = create_cluster_based_data_structure(tissue_specific_junction_file, 20)


# Loop through clusters to save
for cluster_id in cluster_struct.keys():
	geno_0_indis = np.asarray([cluster_struct[cluster_id]['geno_0_indi']])
	geno_1_indis = np.asarray([cluster_struct[cluster_id]['geno_1_indi']])
	geno_2_indis = np.asarray([cluster_struct[cluster_id]['geno_2_indi']])

	jxn_names = cluster_jxn_data_structure[cluster_id]['jxn_names']
	jxn_matrix = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
	ordered_individuals = np.asarray(cluster_jxn_data_structure[cluster_id]['samples'])

	# Save results for inliers
	inlier_output_file = visualize_cluster_distribution_dir + tissue_name + '_' + cluster_id + '_counts.txt'
	save_jxns_to_output_file(geno_0_indis, geno_1_indis, geno_2_indis, jxn_names, jxn_matrix, ordered_individuals, inlier_output_file)