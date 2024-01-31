#!/usr/bin/python
import sys
import gc
import os
import re
import pickle
import collections
import multiprocessing
import numpy as np
from optparse import OptionParser, OptionGroup
from time import time

from . import gqv_utils
from . import gqv_dbg_utils
from . import gqv_bam_utils
from . import gqv_glm
from . import gqv_software_management as SOFTMAN
from . import gqv_vartool
from . import gqv_gff_utils 
from . import gqv_phaser


__author__  = 'Giovanni Quinones Valdez'


def process_read_cluster(RC_info):
	t0 = time()

	(GeneName, CHROM, RCG_Start, RCG_End, STRAND), RCB_list = RC_info

	genome_faidx = gqv_utils.read_genome_fasta(options.Genome)

	RCG_reads, RCB_reads = gqv_bam_utils.get_RC_partial_reads(options, RC_info, genome_faidx)

	VAR_LIST_RC = collections.defaultdict(dict)
	REF_LIST_RC = collections.defaultdict(dict)

	n_RCBs  = len(RCB_reads)
	n_Reads = sum(len(attr[1]) for rcb, attr in RCB_reads.items()) 

	SOFTMAN.print_time_stamp(f"\tRCG: {GeneName} {CHROM}:{RCG_Start}-{RCG_End} RCBs = {n_RCBs} pReads = {n_Reads}") 

	for RCB_i, (RCB, RCB_attr) in enumerate(RCB_reads.items()):
	
		(RCB_Start, RCB_End) = RCB
		(REFPOS, PartialReads) = RCB_attr

		# SOFTMAN.print_time_stamp(f"\t\tRCB: {RCB_i} {RCB_Start}-{RCB_End} cov = {len(PartialReads) - 1}")
		### Construct de Bruijn graph
		db = gqv_dbg_utils.construct_de_Bruijn_graph(PartialReads, options) 
		
		### Process de Bruijn graph
		db.name = f"RC:{RCB_Start}-{RCB_End}" 

		xdb, SetOfReadIndexes = gqv_utils.Process_de_Bruijn_Graph(PartialReads, options, db)  
		del db

		REF_INDEX = PartialReads[0]['Index'] 

		# gqv_dbg_utils.draw_de_Bruijn_graph(xdb, options.kmer, f'{RCB_i}.dbg', Ref_Index = REF_INDEX)

		### DFT algorithm for src-snk pairs compressed
		Src_Snk_pairs = gqv_dbg_utils.find_Source_Target_pairs3(xdb, REF_INDEX, SetOfReadIndexes, options) 

		del xdb

		### Find variants for every read index
		SetOfVariants = gqv_utils.assign_vars_to_reads(Src_Snk_pairs, REF_INDEX, SetOfReadIndexes, 
														PartialReads, options.kmer, REFPOS)  

		SOFTMAN.rm_nested_dict(Src_Snk_pairs) 
		
		### Process, filter and modify variants
		SetOfVariants = gqv_vartool.find_equivalent_indels(SetOfVariants, REFPOS) 

		### Obtian all overlapping read indexes for every variant
		SetOfVariants, REF_blocks = gqv_utils.overlap_vars_and_reads(SetOfVariants, SetOfReadIndexes, 
														   PartialReads, REFPOS)  

		# SOFTMAN.rm_nested_dict(PartialReads)
		SOFTMAN.rm_nested_list(SetOfReadIndexes) 
		
		### Feature collection for every variant
		VAR_LIST = gqv_utils.get_variant_attr(SetOfVariants, REFPOS, GeneName)
		REF_LIST = gqv_utils.get_reference_attr(REF_blocks, REFPOS, GeneName)

		SOFTMAN.merge_copy(VAR_LIST_RC, VAR_LIST)  
		SOFTMAN.merge_copy(VAR_LIST_RC, REF_LIST) 


	RCG_reads_list = list(RCG_reads.keys())

	t7 = time()

	VAR_LIST_RC, hap_clusters = gqv_phaser.main(VAR_LIST_RC, RCG_reads_list, clf_dict, options)

	t8 = time()

	gene_pickle_file = f'{options.gff_file}/tx_structure.{GeneName}.pickle'
	

	if options.run_mode == "Training":
		mi_outlines = []
	elif not os.path.exists(gene_pickle_file):
		mi_outlines = []
	elif not hap_clusters:
		mi_outlines = []
	else:
		tb = open(gene_pickle_file, 'rb')
		TX_structure = pickle.load(tb)

		mi_outlines = gqv_gff_utils.mi_parse_variants(VAR_LIST_RC, RCG_reads, TX_structure, 
														RC_info, hap_clusters, ami_polyfits, ami_normalfits)
		del TX_structure
		tb.close()
		
	t9 = time()

	SOFTMAN.print_time_stamp(f"\tRCG:  {CHROM}:{RCG_Start}-{RCG_End} tx_str = {t9-t8:.2e} phaser = {t8-t7:.2e} main = {t7-t0:.2e} done")
	SOFTMAN.merge_copy(VAR_LIST_RC, REF_LIST_RC)

	return VAR_LIST_RC, mi_outlines




def parse_chroms(options):
	region       = options.search_region
	genome_faidx = gqv_utils.read_genome_fasta(options.Genome)
	bam_chroms   = gqv_bam_utils.get_chroms_from_bam(options.input_bam)

	common_chroms = set(genome_faidx.index) & set(bam_chroms)

	Region_List = []

	r = re.match("([\d|\w]+):(\d+)\-(\d+)", region)

	if not common_chroms:
		# Unmatched assemblies
		raise ValueError("Invalid genome reference\n")
	elif region == "All":
		# Default: all chromosome
		Region_List = [(chrom, None, None) for chrom in common_chroms]
	elif bool(r):
		# Region is an interval
		chrom, fetch_start, fetch_end = r.groups()
		if chrom not in common_chroms:
			raise ValueError(f"Region is not valid {region}\n")
		Region_List = [(chrom, int(fetch_start), int(fetch_end))]
	else:
		# Region is a chromosome
		if region in common_chroms:
			Region_List = [(region, None, None)]
		else:
			raise ValueError(f"Region is not valid {region}\n")

	Region_List.sort()
	
	return Region_List



def main():

	SOFTMAN.print_time_stamp("COMMAND = '{}'\n".format(' '.join(sys.argv)))

	description = """isoLASER"""

	usage = """\n\tpython {} -b <file.bam> -g <genome.fa> -t <tx_str_dir> -o <prefix>""".format(sys.argv[0])

	parser = OptionParser(usage = usage, description = description)
	
	parser.add_option("-b", "--input-bam",   
		dest = "input_bam", 
		help = "[Required] Input bam file, sorted and indexed")
	parser.add_option("-o", "--output-prefix",  
		dest = "output_prefix", 
		help = "[Required] Prefix of the output files")
	parser.add_option("-f", "--genome-file", 
		dest = "Genome", 
		help = "[Required] Reference Genome (Fasta format)")
	parser.add_option("-t", "--gff-file", 
		dest = "gff_file", 
		help = "[Required] Transcript annotation (GFF format)")
	parser.add_option("-s", "--sample-name",  
		dest = "sample_name", 
		help = "[Optional] sample name", default = None)


	filtering_options = OptionGroup(parser, "Filtering Options")

	filtering_options.add_option("--AB", 
		dest = "minRatioVar", 
		type = "float",  
		default = 0.05, 
		help = "Minimum allelic ratio for the variant allele. [Default: 0.05]")
	filtering_options.add_option("--AC", 
		dest = "minCountVar", 
		type = "int",    
		default = 3,    
		help = "Minimum read depth supporting the variant allele. [Default: 2]")
	filtering_options.add_option("--DP", 
		dest = "minCoverage", 
		type = "int",    
		default = 10,   
		help = "Minimum read depth at the position of the variant. [Default: 10]")
	
	parser.add_option_group(filtering_options)
	

	run_mode_options = OptionGroup(parser, "Run Mode Options")

	run_mode_options.add_option("--run_mode", 
		dest = "run_mode", 
		default = 'Variant_Caller', 
		help = "Select <Variant_Caller> or <Training> mode.[Default: Variant_Caller]")
	run_mode_options.add_option("--platform", 
		dest = "platform", 
		default = "PacBio", 
		help = "Select from (PacBio, Nanopore). [Default: PacBio]")
	
	parser.add_option_group(run_mode_options)



	advanced_options = OptionGroup(parser, "Advanced Options")

	advanced_options.add_option("-c", "--region", 
		dest = "search_region",     
		default = "All", 
		help = "Limit search to this region (CHROM:start-end or CHROM). [Default: All]")
	advanced_options.add_option("--gene", 
		dest = "gene_id",     
		default = None, 
		help = "Limit search to this gene id. [Default: None]")
	advanced_options.add_option("-n", "--nodes", 
		dest = "nodes", 
		type = "int", 
		default = 64, 
		help = "Number of threads for parallel execution. [Default = 64]")
	advanced_options.add_option("--min-map_quality", 
		dest = "minMapQ", 
		type = "int", 
		default = 20,
		help = "Minimum mapping quality to keep a read. [Default = 20]")
	advanced_options.add_option("--include-continuous", 
		dest = "include_continuous", 
		action = "store_true",
		default = False,
		help = "Include reads with no splice junctions (`N` cigar). [Default = False]")
	advanced_options.add_option("--kmer-size", 
		dest = "kmer", 
		type = "int", 
		default = 21,
		help = "k-mer size for de-Bruijn graph assembly. [Default: 15]")
	advanced_options.add_option("--maxSeqErrorRate", 
		dest = "maxSeqErrorRate", 
		type = "float", 
		default = 1e-1,
		help = "Maximum estimate of sequencing error rate. [Default: 0.1]")
	advanced_options.add_option("--Ploidy", 
		dest = "ploidy", 
		type = "int", 
		default = 2,
		help = "Maximum ploidy to be considered. [Default: 2]")


	parser.add_option_group(advanced_options)

	global options

	(options, args) = parser.parse_args()

	sys.setrecursionlimit(int(1e6))

	if not os.path.exists(options.input_bam):
		parser.error('Bam file does not exist')

	if not os.path.exists(options.Genome):
		parser.error('Genome file does not exist')

	if not options.output_prefix:   
		parser.error('Output file prefix not given')
	
	if options.kmer > 25 or options.kmer < 15:
		SOFTMAN.print_time_stamp("WARNING: We recomend a k-mer size of 15-25")


	PACKAGE_PATH = os.path.realpath(os.path.dirname(__file__))

	global clf_dict, ami_polyfits, ami_normalfits

	if options.platform == "PacBio":
		with open(f"{PACKAGE_PATH}/glm/gm12878_sequel2.MERGED.feature_matrix.tab.proc_1.csv.f1.glm.pickle", 'rb') as fid:
			clf_dict = pickle.load(fid)
	elif options.platform == "Nanopore":
		with open(f"{PACKAGE_PATH}/glm/merged.GTEX.MERGED.ALLCHR.feature_matrix.tab.proc_1.csv.f1.glm.pickle", 'rb') as fid:
			clf_dict = pickle.load(fid)
	else:
		parser.error("Platform not valid {}. Select from [PacBio, Nanopore]\n".format(options.platform))	

	
	with open(f'{PACKAGE_PATH}/glm/ami_null.polyfits.ar_psi_bins.pickle', 'rb') as f:
		ami_polyfits = pickle.load(f)
	with open(f'{PACKAGE_PATH}/glm/mi_testing.training_summary.normalfits.pickle', 'rb') as f:
		ami_normalfits = pickle.load(f)

	if options.run_mode not in ("Variant_Caller", "Training"):
		parser.error("Run mode {} not valid\n".format(options.run_mode))


	#--------- Run info ----------
	n = min(options.nodes, multiprocessing.cpu_count())
	SOFTMAN.print_time_stamp(f"Running isoLASER \n")
	SOFTMAN.print_time_stamp(f"Running {n} child processes\n")


	Search_Regions = parse_chroms(options)

	
	# Region search (usually chromosome)
	write_header = True 

	for (CHROM, fetch_start, fetch_end) in Search_Regions:

		CHROM_VARS_FEAT   = collections.defaultdict(dict)
		ALL_VARS_MUTINFO  = []

		GENE_LIST = gqv_gff_utils.get_genes_from_gtf(options.gff_file, CHROM, fetch_start, fetch_end)

		if options.gene_id is not None and options.gene_id not in GENE_LIST:
			continue 

		RCG = gqv_bam_utils.find_read_clusters(options, CHROM, fetch_start, fetch_end, GENE_LIST)

		if RCG:
			SOFTMAN.print_time_stamp(f"Contig {CHROM} : {len(RCG)} gene clusters found.")

		if options.gene_id is not None:
			RCG = [rcg for rcg in RCG if rcg[0][0] == options.gene_id]

		pool = multiprocessing.Pool(n) 
		pool_output = pool.imap_unordered(process_read_cluster, RCG)

		for (RC_VAR_LIST, RC_MI_LIST) in pool_output: 
			SOFTMAN.merge_copy(CHROM_VARS_FEAT, RC_VAR_LIST)
			ALL_VARS_MUTINFO.extend(RC_MI_LIST)

		pool.terminate()
			
			

		if options.run_mode == "Variant_Caller":
			SM = gqv_bam_utils.get_sample_name(options.input_bam, options.sample_name)

			gqv_utils.write_vcf_file(CHROM_VARS_FEAT, Search_Regions, options.Genome,
										options.output_prefix + ".gvcf", SM, options,
										write_header = write_header)

			gqv_utils.write_mutinfo_file(ALL_VARS_MUTINFO, 
										options.output_prefix + ".mi_summary.tab",
										write_header = write_header)

		elif options.run_mode == "Training":
			gqv_utils.write_feat_file(CHROM_VARS_FEAT, 
										options.output_prefix + ".feature_matrix.tab",
										write_header = write_header)

		write_header = False 


	SOFTMAN.print_time_stamp("Job Done")
	SOFTMAN.print_time_stamp(f"Output prefix: {options.output_prefix}")



if __name__ == "__main__":
	main()





