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


def fetch_gtf_records(gtf, gene_id, chrom, POS, padding = 1500):
	REGION = [POS - padding, POS + padding]

	gene_convert_id = {} 

	for record in gtf.fetch(chrom, POS - padding, POS + padding, parser = pysam.asGTF()):
		if record.feature == "transcript" and record.gene_id.startswith(gene_id):
			if record.source != "TALON":
				if record.start < REGION[0]:
					REGION[0] = record.start
				if record.end > REGION[1]:
					REGION[1] = record.end
		elif record.feature == "gene":
			gene_convert_id[record.gene_id] = record.gene_name

	REGION = f'{chrom}:{REGION[0]}-{REGION[1]}'

	try:
		gene_name = gene_convert_id[gene_id]
	except:
		gene_name = gene_id

	return REGION, gene_name


def read_target_gene_file(gene_file):
	target_genes = set()
	with open(gene_file) as f:
		for line in f:
			gene_id, gene_name = line.strip().split('\t')
			gene_id = gene_id.split('.')[0] # remove gene id version
			target_genes.add(gene_id)

	return target_genes



if __name__ == "__main__":

    logging.basicConfig(level = logging.INFO)

    parser = optparse.OptionParser()

    parser.add_option('--mi', dest = 'mi_file', \
        help = 'mi_file')
    parser.add_option('-o', dest = 'output_file', \
        help = 'output file', default = None)
    parser.add_option('-g', dest = 'GTF', \
        help = 'gtf file', default = None)
    parser.add_option('--genes', dest = 'gene_file', \
        help = 'target genes', default = None)

    options, args = parser.parse_args()

    if options.output_file is not None:
        taget_genes = read_target_gene_file(options.gene_file)


    ### read gtf
    gtf_handle = pysam.TabixFile(options.GTF)

    ## output files
    bedf = open(options.output_file + '.bed', 'w')
    genf = open(options.output_file + '.gen', 'w')

    var_info = defaultdict(set)


    df = pd.read_csv(options.mi_file, sep = '\t')


    df = df[(df.label == "cis")]


    for row in df.itertuples():
        delta_psi = float(row._6)

        if abs(delta_psi) < 0.05 :
            continue

        gene_id = row.Gene_name
        gene_id_versionless = gene_id.split('.')[0]

        if options.gene_file is not None and gene_id_versionless not in taget_genes:
            continue

        if row.COV >= 20:

            chrom, pos, alleles = row.variant_qual.split(':')
            ref, alt = alleles.split('>')

            if len(ref) > 1 or len(alt) > 1:
                continue #indel 

            var_info[gene_id].add((row.variant_qual, abs(delta_psi)))

            exonic_part = re.match(r'(\w+):\[(\d+),(\d+)\)/(\S+):(\S+)', row.exonic_part)
            bed_line = [exonic_part.group(1), exonic_part.group(2), exonic_part.group(3), gene_id, '.']
            bedf.write('\t'.join(map(str, bed_line)) + '\n')


    for gene, var_set in var_info.items():
           
        (best_var, best_dpsi) = max(var_set, key = lambda x: x[1])
        chrom, pos, alleles = best_var.split(':')

        gene_interval, gene_name = fetch_gtf_records(gtf_handle, gene, chrom, int(pos))

        gene_line = [gene, gene_name, best_var, gene_interval]
        genf.write('\t'.join(gene_line) + '\n')

