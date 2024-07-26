#!/usr/bin/python
import sys
import pysam
import HTSeq 
import numpy as np
import collections
import pickle
from optparse import OptionParser, OptionGroup
import os
import os.path



def reads_by_gene(bam_file):

	transcriptome_map = {}

	bh = pysam.AlignmentFile(bam_file)

	for read in bh:
		if not read.is_secondary and not read.is_supplementary:
			transcriptome_map[read.query_name] = read.reference_name

	return transcriptome_map 



def read_by_assignment(assign_file):
	
	assignment_map = {}

	with open(assign_file, 'r') as fh:
		for line in fh:
			line = line.strip().split('\t')
			assignment_map[line[0]] = line[1]

	return assignment_map



def main():
	description = "isoLASER - filter and annotate bam file"

	usage = ""

	parser = OptionParser(usage = usage, description = description)
	
	parser.add_option("-b", "--input-bam",   
		dest = "input_bam", help = "[required]")

	parser.add_option("-o", "--output-bam",  
		dest = "output_bam", help = "[required]")
	parser.add_option("-g", "--gtf_file", 
		dest = "gtf_file", help = "[required]")
	
	annotation_group = OptionGroup(parser, "Annotation Options")

	annotation_group.add_option("-t", "--tx-bam",   
		dest = "transcriptome_bam")
	annotation_group.add_option("-a", "--assignment-file",   
		dest = "assignment_file")

	parser.add_option_group(annotation_group)
	
	parser.add_option("--padding", 
		dest = "padding", type = int, default = 10)

	(options, args) = parser.parse_args()

	if options.transcriptome_bam is not None:
		sys.stderr.write(f"Reading transcriptome bam file: {options.transcriptome_bam}\n")
		transcriptome_map = reads_by_gene(options.transcriptome_bam)
	elif options.assignment_file is not None:
		sys.stderr.write(f"Reading assignment file: {options.assignment_file}\n")
		transcriptome_map = read_by_assignment(options.assignment_file)
	else:
		raise ValueError(f"Transcriptome information not provided\n")

	bh = pysam.AlignmentFile(options.input_bam)
	oh = pysam.AlignmentFile(options.output_bam, 'wb', template = bh)


	r_stack = set()	
	r_counter = collections.Counter()
	g_counter = 0 

	for record in HTSeq.GFF_Reader(options.gtf_file):
		
		if record.type != "gene":
			continue 

		gene_name = record.attr["gene_id"] + '_' + record.attr["gene_name"]
		gene_name = gene_name.replace('"', '')
		g_counter += 1

		if g_counter % 2500 == 0:
			sys.stderr.write(f'\tgenes processed: {g_counter:,}. reads processed: {r_counter["total"]:,}\n')
			sys.stderr.flush()

		for read in bh.fetch(record.iv.chrom, record.iv.start, record.iv.end):
			r_counter['total'] += 1

			r_name = read.query_name 	

			### filter if read is secondary
			if read.is_secondary or read.is_supplementary or r_name in r_stack:
				continue 

			r_counter['basic_filter'] += 1

			### filter if read is trans-gene
			if read.positions[0] < (record.iv.start - options.padding) or read.positions[-1] > (record.iv.end + options.padding):
				continue 
			
			r_counter['trans_gene'] += 1

			if r_name not in transcriptome_map:
				continue

			r_counter['annotated'] += 1

			read.set_tag('ZG', gene_name)

			if transcriptome_map[r_name]:
				read.set_tag('ZT', transcriptome_map[r_name])
			else:
				read.set_tag('ZT', 'nan')

			oh.write(read)

			r_stack.add(r_name)

			

	bh.close()
	oh.close() 	

	for filter_value in ['total', 'basic_filter', 'trans_gene', 'annotated']:
		sys.stderr.write(f'\t{filter_value}: {r_counter[filter_value]}\n')

if __name__ == "__main__":
	main()