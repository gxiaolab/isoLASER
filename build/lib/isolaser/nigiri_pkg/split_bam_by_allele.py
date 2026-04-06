#!/usr/bin/python
import pysam
import pandas as pd
import sys
import re
import os
import logging
import optparse
from collections import Counter, defaultdict
import subprocess 
from scipy.stats import chisquare, norm, chi2_contingency

from . import gqv_bam_utils 


def process_hla_sam(hla_sam):
	d = {}
	bh = pysam.AlignmentFile(hla_sam)

	for record in bh:
		if record.flag == 0:
			try: 
				allele_group = record.reference_name.split(':')[1].split('_')[1]
			except: 
				allele_group = 'na'
			read_name = record.query_name 
			d[read_name] = allele_group
	return d 


def main():

	logging.basicConfig(level = logging.INFO)

	parser = optparse.OptionParser()

	parser.add_option('-b', dest = 'bam', \
		help = 'bam file (sorted and indexed)')
	parser.add_option('-v', dest = 'var', \
		help = 'variant string `CHROM:POS0-POS1:REF>ALT`')
	parser.add_option('--hla', dest = 'hla_sam', \
		help = 'hla aligned reads', default = None)
	parser.add_option('-o', dest = 'output_file', \
		help = 'output file', default = None)
	parser.add_option('-s', dest = 'sample_name', \
		help = 'sample name', default = 'sample_1')
	parser.add_option('-t', dest = 'tmpdir', \
		help = 'output file', default = 'tmp')

	options, args = parser.parse_args()

	TMPDIR = options.tmpdir

	tx_allele_counter = defaultdict(lambda : defaultdict(lambda : 0))

	### check hla reads
	if options.hla_sam is not None:
		reads_hla_dict = process_hla_sam(options.hla_sam)
	else:
		reads_hla_dict = {} 


	# process var
	chrom, POS0, REF, ALT = re.split(':|>', options.var)

	POS = int(POS0) 
	END = int(POS0)

	logging.info(f"splitting by locus: {chrom}:{POS}, alleles = [{REF}, {ALT}]")


	bam_handle = pysam.Samfile(options.bam, 'rb')

	ref_bam = '{}/allele={}.bam'.format(TMPDIR, REF)
	alt_bam = '{}/allele={}.bam'.format(TMPDIR, ALT)

	logging.info(f"writting to {ref_bam}\n{alt_bam}")

	bam_allele_1 = pysam.Samfile(ref_bam, 'wb', bam_handle)
	bam_allele_2 = pysam.Samfile(alt_bam, 'wb', bam_handle)

	read_counter = Counter()

	# HLA typing
	HLA_TYPING = defaultdict(lambda : defaultdict(lambda : 0))

	

	for column in bam_handle.pileup(chrom, POS - 2, POS + 1, truncate = True):
		if column.pos == POS:
			for r in column.pileups:
				if r.is_del or r.is_refskip:
					continue
				
				base  = r.alignment.query_sequence[r.query_position]
				tx_id = r.alignment.get_tag('ZT')

				try:
					hla_contig = reads_hla_dict[r.alignment.query_name]
				except:
					hla_contig = 'allele'

				if base == REF:
					is_bad, lab = gqv_bam_utils.filter_reads(r.alignment)
					read_counter[('ref', lab)] += 1
					if not is_bad:
						tx_allele_counter[REF][tx_id] += 1
						bam_allele_1.write(r.alignment)
						
						HLA_TYPING[REF][hla_contig] += 1

				elif base == ALT:
					is_bad, lab = is_bad_read(r.alignment)
					read_counter[('alt', lab)] += 1
					if not is_bad:
						tx_allele_counter[ALT][tx_id] += 1
						bam_allele_2.write(r.alignment)

						HLA_TYPING[ALT][hla_contig] += 1

	bam_allele_1.close()
	bam_allele_2.close()

	logging.info(f"sorting and indexing bams")
	logging.info(f"\ttotal reference allele {REF} reads: {read_counter}")
	
	tx_allele_counter = pd.DataFrame(tx_allele_counter)
	tx_allele_counter.fillna(0, inplace = True)
	tx_allele_counter = tx_allele_counter.astype(int)

	real_counts = [list(tx_allele_counter[REF]), list(tx_allele_counter[ALT])]
	chi2, p_value, a, b = chi2_contingency(real_counts)



	logging.info(f"\ttx alleles counter:\n{tx_allele_counter}")
	logging.info(f"\tsm {options.sample_name} asts p-value= {p_value:.3e}")
	#logging.info(f"\ttotal alternative allele {ALT} reads: {read_counter[2]} passed: {read_counter[3]}")
	
	sorted_bam_ref = ref_bam.replace('.bam', '.sorted.bam')
	sorted_bam_alt = alt_bam.replace('.bam', '.sorted.bam')

	pysam.sort("-o", sorted_bam_ref, ref_bam)
	pysam.index(sorted_bam_ref)

	pysam.sort("-o", sorted_bam_alt, alt_bam)
	pysam.index(sorted_bam_alt)

	os.remove(ref_bam)
	os.remove(alt_bam) 

	### write fofn
	if options.output_file is not None:
		
		HLA_REF = max(HLA_TYPING[REF], key = HLA_TYPING[REF].get)
		HLA_ALT = max(HLA_TYPING[ALT], key = HLA_TYPING[ALT].get)

		with open(options.output_file, 'w') as f:
			f.write(f'id1\t{os.path.abspath(sorted_bam_ref)}\t{REF}-{HLA_REF}\n')
			f.write(f'id2\t{os.path.abspath(sorted_bam_alt)}\t{ALT}-{HLA_ALT}\n')


	logging.info(f"all done")



if __name__ == "__main__":
	main() 
