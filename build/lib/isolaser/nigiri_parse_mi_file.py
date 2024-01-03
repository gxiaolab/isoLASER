#!/usr/bin/python
import pandas as pd
import sys
import re
import os
import logging
import optparse
from collections import Counter, defaultdict
import subprocess 
import pysam

from . import gqv_gff_utils



def main():

	logging.basicConfig(level = logging.INFO)

	parser = optparse.OptionParser()

	parser.add_option('--mi', dest = 'mi_file', \
		help = 'mi_file')
	parser.add_option('-o', dest = 'output_file', \
		help = 'output file', default = None)
	parser.add_option('-t', dest = 'transcript_db', \
		help = 'gtf file', default = None)
	
	parser.add_option('--genes', dest = 'gene_file', \
		help = 'target genes', default = None)

	options, args = parser.parse_args()


	### read genes
	gene_lst = gqv_gff_utils.get_genes_from_gtf(options.transcript_db)

	## output files
	bedf = open(options.output_file + '.cis_events.bed', 'w')
	genf = open(options.output_file + '.cis_genes.tsv', 'w')

	var_info  = defaultdict(set)
	cis_genes = defaultdict(set)

	df = pd.read_csv(options.mi_file, sep = '\t', escapechar='#')

	df = df[(df.label == "cis")]

	for row in df.itertuples():

		if row.Delta_PSI < 0.05:
			continue
		if (row.end - row.start) < 6:
			continue 

		gene_id = row.Gene_name

		exonic_part = (row.chrom, row.start, row.end)

		if row.Linkage_Type == "HAP-ExonicPart":

			cis_genes[gene_id].add(exonic_part)
				
		if row.Linkage_Type == "VAR" and ">" in row.PG:

			if row.Cov >= 20:

				chrom, pos, alleles = row.PG.split(':')
				ref, alt = alleles.split('>')

				if len(ref) > 1 or len(alt) > 1:
					continue 

				if row.AR >= 0.1 and row.AR <= 0.9 and row.PSI >= 0.1 and row.PSI <= 0.9:
					variant = row.PG
					var_info[gene_id].add((variant, exonic_part, row.Delta_PSI))

	logging.info(f"writing to:\n\t- {options.output_file}.cis_events.bed\n\t- {options.output_file}.cis_genes.tsv")
	logging.info(f"found {len(set(var_info) & set(cis_genes))}  genes with cis events")

	for gene, var_set in var_info.items():

		if gene not in cis_genes or not var_set:
			continue
		
		var_set = list(var_set)

		for (variant, exonic_part, dpsi) in var_set[:]:
			
			if exonic_part not in cis_genes[gene]:
				var_set.remove((variant, exonic_part, dpsi))
			else:
				
				bed_line = list(exonic_part) + [gene_id, dpsi]
				bedf.write('\t'.join(map(str, bed_line)) + '\n')


		(best_var, best_xp, best_dpsi) = max(var_set, key = lambda x: x[2])
		chrom, pos, alleles = best_var.split(':')

		gene_interval = "{0}:{1}-{2}".format(*gene_lst[gene_id])

		gene_line = [gene, best_var, gene_interval, dpsi]
		genf.write('\t'.join(map(str, gene_line)) + '\n')

	bedf.close()
	genf.close()
	logging.info(f"done")


if __name__ == "__main__":
	main() 