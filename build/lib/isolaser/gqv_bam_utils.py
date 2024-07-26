#!/usr/bin/python
import pysam
from collections import defaultdict
from HTSeq import GenomicInterval
import numpy as np
import random
import sys


def filter_reads(pysam_read, P_minMapQual, include_continuous):
	tags = dict(pysam_read.tags)

	if pysam_read.mapping_quality < P_minMapQual:
		return 'mapq', True
	elif len(pysam_read.query_alignment_sequence) < 100:
		return 'rlen', True
	elif len(pysam_read.query_alignment_sequence)/len(pysam_read.query_sequence) < 0.67:
		return 'sclip', True
	elif "N" in pysam_read.query_alignment_sequence:
		return 'N', True
	elif pysam_read.is_unmapped:
		return 'unmap', True
	elif pysam_read.is_secondary:
		return 'sec', True
	elif pysam_read.is_duplicate:
		return 'dup', True
	elif pysam_read.is_qcfail:
		return 'qcfail', True
	elif "NH" in tags and tags["NH"] > 1:
		return 'mmap', True
	elif pysam_read.is_supplementary:
		return 'supp', True
	elif pysam_read.get_tag('ZT') == 'nan':
		return 'zt_nan', True
	# elif pysam_read.get_tag('ZX') in ('Intergenic', 'Antisense', 'Genomic', 'ISM'):
	# 	return 'zx_bad', True
	elif "N" not in pysam_read.cigarstring and not include_continuous:
		return 'no_jxn', True
	else:
		return 'good', False



def process_CIGAR(genome_pos, cigar, more_blocks = None):
	genome_block = 0
	genome_block_lengths = [] 

	read_block = 0
	read_block_lengths = [] 
	read_pos = 0

	matched_blocks = []

	for (ctype, length) in cigar:
		if ctype in (0, 7, 8):   # M,=,X	
			mb = [(read_pos + read_block, 
				   read_pos + read_block + length), 
				  (genome_pos + genome_block, 
				   genome_pos + genome_block + length)]
			matched_blocks.append(mb)
			genome_block += length
			read_block   += length
		elif ctype == 1: #I
			read_block   += length
		elif ctype == 2: #D
			genome_block += length
		elif ctype == 3: #N
			read_block_lengths.append((read_pos, read_pos + read_block))
			read_pos += read_block
			read_block = 0

			genome_block_lengths.append((genome_pos, genome_pos + genome_block))
			genome_pos += (genome_block + length)
			genome_block = 0
	
	genome_block_lengths.append((genome_pos, genome_pos + genome_block))
	read_block_lengths.append((read_pos, read_pos + read_block))

	return read_block_lengths, genome_block_lengths, matched_blocks




def get_chroms_from_bam(input_bam):
	bam_handle = pysam.Samfile(input_bam, 'rb')
	bam_chroms = bam_handle.references
	bam_handle.close()
	return bam_chroms



def get_sample_name(input_bam, sample_name):
	if sample_name is not None:
		return sample_name
	else:
		bam_handle = pysam.Samfile(input_bam, 'rb')	
		try:
			return  bam_handle.header['RG'][0]['SM']
		except KeyError:
			return  input_bam



def find_read_clusters(options, CHROM, start, end, GENE_LIST):
	
	read_block_coords = defaultdict(lambda : defaultdict(lambda : 0))

	bam_handle = pysam.Samfile(options.input_bam, 'rb')

	log_reads_filtered = defaultdict(int)

	for read in bam_handle.fetch(CHROM, start = start, end = end):
		Label, no_pass = filter_reads(read, options.minMapQ, options.include_continuous)
		log_reads_filtered[Label] += 1
		if no_pass:
			continue

		r_blocks, g_blocks, m_blocks = process_CIGAR(read.pos, read.cigar)
		gene_id = read.get_tag('ZG')
		
		for sB, eB in g_blocks:
			read_block_coords[gene_id][sB] += 1
			read_block_coords[gene_id][eB] -= 1
			
	bam_handle.close()

	for skip_read_label, skip_read_count in log_reads_filtered.items():
		sys.stderr.write(f"[Warning] Reads filtered by {skip_read_label} = {skip_read_count}\n")

	Read_Clusters_Genes = defaultdict(list)

	for gene_id, coord_dict in read_block_coords.items():

		if gene_id not in GENE_LIST:
			nr = sum([x for x in coord_dict.values() if x > 0])
			s = sorted(coord_dict.keys())
			ss, ee = s[0], s[-1] 
			sys.stderr.write(f"[Warning] gene id {gene_id} not found in {options.gff_file} reads = {nr} coord = {ss}-{ee}\n")
			gene_key = (gene_id, CHROM, ss, ee, '.')
		else:
			(gene_iv_chrom, gene_iv_start, gene_iv_end, gene_iv_strand) = GENE_LIST[gene_id]
			gene_key = (gene_id, gene_iv_chrom, gene_iv_start, gene_iv_end, gene_iv_strand)
		
		Cum_Cov_block = 0
		Prv_Cov_block = 0

		coords_cov = sorted(coord_dict.items())
		
		if not coords_cov:
			continue 
		
		prev_coord = coords_cov[0][0]
		_max_cov = 0 
		
		for coord, cov in coords_cov: 
			Cum_Cov_block += cov			
			log10FC = np.log10(Cum_Cov_block + 1) - np.log10(Prv_Cov_block + 1)
			log2FC = np.log2(Cum_Cov_block + 1) - np.log2(Prv_Cov_block + 1)

			if abs(log2FC) >= 1.0 or (Cum_Cov_block * Prv_Cov_block == 0) :	
				Read_Clusters_Genes[gene_key].append((prev_coord, coord))
				prev_coord = coord
			
			Prv_Cov_block = Cum_Cov_block
			if Cum_Cov_block > _max_cov:
				_max_cov = Cum_Cov_block

		Read_Clusters_Genes[gene_key].append((prev_coord, coord))

	return sorted(Read_Clusters_Genes.items())



class genome_ref_pos():
	
	def __init__(self, CHROM, RCB_start, RCB_end, STRAND, genome, k_ext):

		RC_ext_seq_full  = "X"*k_ext
		RC_ext_seq_full += str(genome.fetch(CHROM, RCB_start + 1, RCB_end)).upper()
		RC_ext_seq_full += "Y"*k_ext

		self.Seq_adapter = str(genome.fetch(CHROM, RCB_start - k_ext + 1, RCB_end + k_ext)).upper()

		self.Seq_ext = RC_ext_seq_full
		self.gi = GenomicInterval(CHROM, RCB_start - k_ext, RCB_end + k_ext)
		self.strand = STRAND
		self.k_ext = k_ext

	def get_sequence(self, g_start, g_end, adapter = True):
		off = g_start - self.gi.start 
		assert off >= 0 
		if adapter :
			return self.Seq_ext[off : off + (g_end - g_start)]
		else:
			return self.Seq_adapter[off : off + (g_end - g_start)]

	def genome_pos(self, relative_pos):
		return self.gi.start + relative_pos 




def get_RC_partial_reads(options, RC_info, genome):
	(Gene_Id, CHROM, RCG_Start, RCG_End, STRAND), RCB_list = RC_info

	RCG_reads = defaultdict(dict)

	bam_handle = pysam.Samfile(options.input_bam, 'rb')

	for read in bam_handle.fetch(CHROM, max(0, RCG_Start - 1), RCG_End + 1):

		if filter_reads(read, options.minMapQ, options.include_continuous)[1]:
			continue
		if read.get_tag('ZG') != Gene_Id:
			continue 

		r_blocks, g_blocks, m_blocks = process_CIGAR(read.pos, read.cigar)

		assert len(g_blocks) == len(r_blocks)

		RCG_reads[read.query_name]['read_blocks'] = r_blocks
		RCG_reads[read.query_name]['geno_blocks'] = g_blocks
		RCG_reads[read.query_name]['dir'] = read.is_reverse
		RCG_reads[read.query_name]['seq'] = read.query_alignment_sequence
		RCG_reads[read.query_name]['ZT']  = read.get_tag('ZT')
		# RCG_reads[read.query_name]['ZX']  = read.get_tag('ZX')
						
	bam_handle.close()

	k_ext = (options.kmer - 2)
	RCB_reads = defaultdict(list)

	if len(RCG_reads) > 2500:
		_reads = random.sample(RCG_reads.keys(), 2000)
		for _read in list(RCG_reads):
			if _read not in _reads:
				RCG_reads.pop(_read)

	for RC_i in range(len(RCB_list)):

		PartialReads = defaultdict(dict)
		RCB_start, RCB_end = RCB_list[RC_i] ; print("tmp_rcb", RCB_start, RCB_end)

		for read_name, read_attr in RCG_reads.items():
			for g_i, geno_block in enumerate(read_attr['geno_blocks']):
				if geno_block[0] < RCB_end and RCB_start < geno_block[1]: 
					read_name_partial = f'{read_name}:{RC_i}_{g_i}' 
					sub_r_blocks = list(read_attr['read_blocks'][g_i])
					sub_g_blocks = list(read_attr['geno_blocks'][g_i])
					sub_seq      = read_attr['seq'][sub_r_blocks[0] : sub_r_blocks[1]]

					if sub_g_blocks[0] == RCB_start:
						sub_g_blocks[0] -= k_ext 
						sub_r_blocks[0] -= k_ext 
						sub_seq = 'X'*k_ext + sub_seq

					if sub_g_blocks[1] == RCB_end:
						sub_g_blocks[1] += k_ext 
						sub_r_blocks[1] += k_ext 
						sub_seq = sub_seq + 'Y'*k_ext

					PartialReads[read_name_partial]['seq']         = sub_seq
					PartialReads[read_name_partial]['read_blocks'] = sub_r_blocks 
					PartialReads[read_name_partial]['geno_blocks'] = sub_g_blocks

		if not PartialReads or (RCB_end - RCB_start) < (k_ext + 1):
			continue
		if RCB_start <= k_ext:
			continue 

		REFPOS = genome_ref_pos(CHROM, RCB_start, RCB_end, STRAND, genome, k_ext)
						
		PartialReads[0]['seq']         = REFPOS.Seq_ext
		PartialReads[0]['read_blocks'] = (0, REFPOS.gi.end - REFPOS.gi.start)
		PartialReads[0]['geno_blocks'] = (REFPOS.gi.start, REFPOS.gi.end)

		RCB_reads[(RCB_start, RCB_end)] = (REFPOS, PartialReads)
		

	return RCG_reads, RCB_reads











