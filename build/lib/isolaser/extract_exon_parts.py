#!/usr/bin/python
import HTSeq
import pysam
import numpy as np
import pandas as pd
from optparse import OptionParser, OptionGroup
from collections import defaultdict, Counter

import sys
import os
import multiprocessing 
import pickle 
import time 

from . import gqv_software_management as SOFTMAN


def get_gene_list(gff_file):
    GTF_handle = pysam.TabixFile(gff_file, parser = pysam.asGTF())

    GeneList = []
    for r in GTF_handle.fetch():
        if r.feature == 'gene':
            iv = HTSeq.GenomicInterval(r.contig, r.start, r.end, r.strand)
            GeneList.append((r.gene_id, r.gene_name, iv))
    
    return GeneList



def get_transcript_structure(args):  

    (Gene_id, Gene_Name, Gene_iv) = args 
    CHROM = Gene_iv.chrom

    TX_structure_dict  = defaultdict(dict)

    transcript_dict    = defaultdict(list)

    ### Genomic Arrays
    spl_coords = HTSeq.GenomicArrayOfSets(chroms = [CHROM], stranded = True)
   
    GTF_handle = pysam.TabixFile(GFF_FILE, parser = pysam.asGTF())

    for record in GTF_handle.fetch(CHROM, Gene_iv.start, Gene_iv.end):
        
        if record.gene_id != Gene_id:
            continue
            
        if record.feature == 'exon':
            transcript_dict[record.transcript_id].append((record.start, record.end))


    for transcript_id, exon_list in transcript_dict.items():
        exon_list.sort(key = lambda xn : xn[0])

        for j in range(len(exon_list)):

            if j < len(exon_list) - 1:

                i_s, i_e = (exon_list[j][1], exon_list[j + 1][0])
                intron_iv = HTSeq.GenomicInterval(CHROM, i_s, i_e, Gene_iv.strand)
                spl_coords[intron_iv] += ('intron', transcript_id)

            e_s, e_e = (exon_list[j][0], exon_list[j][1])
            exon_iv = HTSeq.GenomicInterval(CHROM, e_s, e_e, Gene_iv.strand)
            spl_coords[exon_iv] += ('exon', transcript_id)

       
    ### Directionality of the transcript

    for iv, transcript_set in spl_coords.steps(): 

        s = defaultdict(set)
        for (iv_type, tx_id) in transcript_set:
            s[iv_type].add(tx_id)

        if len(s) > 1:
            TX_structure_dict[('exonic_part', str(iv))] = s

    return (Gene_id, Gene_Name, TX_structure_dict)



def main():

    description = """isoLASER - Extracting exon parts from a bam file"""

    usage = """\n\t{} -g <annotation.gtf> -o <prefix>""".format(sys.argv[0])

    parser = OptionParser(usage = usage, description = description)

    parser.add_option("-o", "--output_directory",  
        dest = "output_directory", 
        help = "[Required] Output directory")
    parser.add_option("-g", "--gtf-file", 
        dest = "gtf_file", 
        help = "[Required] GTF file")

    (options, args) = parser.parse_args()


    GeneList = get_gene_list(options.gtf_file)

    # Write Gene List to file
    df = pd.DataFrame(GeneList)
    df.to_csv(f'{options.output_directory}/GeneList.tsv', sep = '\t', index = False, header = False)
            
    SOFTMAN.print_time_stamp(f"Total genes found = {len(GeneList)}") 

    ### multiprocessing pool 

    global GFF_FILE 
    GFF_FILE = options.gtf_file


    pool = multiprocessing.Pool(8) 
    pool_output = pool.imap_unordered(get_transcript_structure, GeneList)

    for g_i, (gene_id, gene_name, gene_tx_structure) in enumerate(pool_output): 

        if g_i % 2000 == 0:
            SOFTMAN.print_time_stamp(f"\t[Progress] {g_i} genes processed")

        gene_pickle_file = f'{options.output_directory}/tx_structure.{gene_id}_{gene_name}.pickle'

        if not gene_tx_structure:
            continue

        with open(gene_pickle_file, 'wb') as pt:
            pickle.dump(gene_tx_structure, pt)

        

    pool.terminate() 



if __name__ == '__main__':
    main()