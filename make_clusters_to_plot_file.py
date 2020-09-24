import numpy as np 
import os
import sys
import pdb
import gzip



def get_individuals_with_rna_seq_and_wgs(genotype_file, junction_file):
	genotype_indi = {}
	merged_indi = {}
	f = open(genotype_file.rstrip())
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for ele in data[2:]:
				genotype_indi[ele] = 1
			continue
	f.close()
	f = gzip.open(junction_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			for ele in data:
				indi = ele.split('-')[0] + '-' + ele.split('-')[1]
				if indi in genotype_indi:
					merged_indi[indi] = 1
			continue
	f.close()
	return merged_indi


def get_variant_position_and_chrom(genotype_file):
	f = open(genotype_file.rstrip())
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		chrom_num = data[0]
		position = data[1]
	f.close()
	return position, chrom_num

def parse_genotypes(genotypes):
	numeric_genotypes = []
	for genotype in genotypes:
		if genotype == './.':
			numeric_genotypes.append(0)
		else:
			numeric_genotype = np.sum(np.asarray(genotype.split('/')).astype(float))
			numeric_genotypes.append(numeric_genotype)
	return numeric_genotypes

def seperate_individuals_by_genotype(genotype_file, valid_individuals):
	f = open(genotype_file.rstrip())
	head_count = 0
	genotype_0 = []
	genotype_1 = []
	genotype_2 = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			indis = data[2:]
			continue
		genotypes = data[2:]
		numeric_genotypes = parse_genotypes(genotypes)
		for i, indi in enumerate(indis):
			if indi not in valid_individuals:
				continue
			genotype = numeric_genotypes[i]
			if genotype == 0.0:
				genotype_0.append(indi)
			elif genotype == 1.0:
				genotype_1.append(indi)
			elif genotype == 2.0:
				genotype_2.append(indi)
			else:
				print('no hit')
				pdb.set_trace()
	f.close()
	return np.asarray(genotype_0), np.asarray(genotype_1), np.asarray(genotype_2)


def get_strand(exon_file, ensamble_id):
	f = open(exon_file)
	strand = 'h'
	for line in f:
		line = line.rstrip()
		data = line.split()
		line_strand = data[3]
		line_ensamble_id = data[4].split('.')[0]
		if line_ensamble_id != ensamble_id:
			continue
		if strand == 'h':
			strand = line_strand
		else:
			if strand != line_strand:
				print('assumption error')
				pdb.set_trace()
	f.close()
	if strand == 'h':
		strand = '+'
	return strand

ensamble_id = sys.argv[1].split('.')[0]
cluster_id = sys.argv[2]
genotype_file = sys.argv[3]
junction_file = sys.argv[4]
exon_file = sys.argv[5]
cluster_to_plot_file = sys.argv[6]



individuals = get_individuals_with_rna_seq_and_wgs(genotype_file, junction_file)

variant_position, variant_chrom = get_variant_position_and_chrom(genotype_file)

genotype_0_indi, genotype_1_indi, genotype_2_indi = seperate_individuals_by_genotype(genotype_file, individuals)

cluster_ids = np.asarray([cluster_id])

strand = get_strand(exon_file, ensamble_id)

t = open(cluster_to_plot_file, 'w')
t.write('cluster_id\tensamble_id\tchrom_num\tvariant_position\tgeno_0_indi\tgeno_1_indi\tgeno_2_indi\tstrand\n')

for cluster_id in cluster_ids:
	t.write(cluster_id + '\t' + ensamble_id + '\t' + variant_chrom + '\t' +variant_position + '\t' + ','.join(genotype_0_indi) + '\t' + ','.join(genotype_1_indi) + '\t' + ','.join(genotype_2_indi) + '\t' + strand + '\n')
t.close()