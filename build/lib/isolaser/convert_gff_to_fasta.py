#!/usr/bin/python

import sys
import Bio
from Bio import SeqIO
import HTSeq 
from pyfaidx import Faidx
import pysam
import pandas as pd
from collections import defaultdict
import optparse



def get_gene_list(GTF_file):
    GeneList = defaultdict(dict)
    GTF_handle = pysam.TabixFile(GTF_file, parser = pysam.asGTF())

    for r in GTF_handle.fetch():
    
        if r.feature == 'exon':
            transcript_bin = (r.contig, r.transcript_id, r.strand)
            try:
                GeneList[r.gene_name][transcript_bin][(r.start, r.end)] = r
            except:
                GeneList[r.gene_name][transcript_bin] = {(r.start, r.end): r}
    return GeneList



def main():
    optparser = optparse.OptionParser(usage = "%prog [options]")

    optparser.add_option("-g", dest = "gtf_in", \
        help = "GTF file", default = None)
    optparser.add_option("-f", dest = "fasta", \
        help = "Fasta file", default = None)   
    optparser.add_option("-o", dest = "fasta_output", \
        help = "Output fasta", default = None)

    (options, args) = optparser.parse_args(sys.argv)

    out_dna_fasta = open(options.fasta_output, 'w')

    FastaFai = Faidx(options.fasta)

    Gene_Dict = get_gene_list(options.gtf_in)

    transcript_count = 0

    for Gene_Name, transcripts_dict in Gene_Dict.items():
        for (chrom, transcript_id, strand), gtf_lines in transcripts_dict.items():

            transcript_count += 1
            if transcript_count % 5000 == 0:
                sys.stderr.write(f'\ttranscripts processed: {transcript_count:,}\n')
                sys.stderr.flush()

            transcript_seq = ""

            for coord, gtf_line in sorted(gtf_lines.items()):
                seq = FastaFai.fetch(chrom, gtf_line.start + 1, gtf_line.end)
                transcript_seq += str(seq) 

            seq = Bio.Seq.Seq(transcript_seq)
            record = Bio.SeqRecord.SeqRecord(seq, id = transcript_id, description="")

            SeqIO.write(record, out_dna_fasta, "fasta")




if __name__ == "__main__":
    main()
           
        
