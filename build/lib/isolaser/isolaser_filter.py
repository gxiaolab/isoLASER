#!/usr/bin/python
import pandas as pd
import numpy as np
import sys
import os
from optparse import OptionParser, OptionGroup
import logging

def main():

    logging.basicConfig(level = logging.INFO)

    description = """isoLASER - filter mi summary file"""

    usage = """\n\tpython {} -m <mi_summary.tsv>  -o <mi_summary.filtered.tsv>""".format(sys.argv[0])

    parser = OptionParser(usage = usage, description = description)

    parser.add_option("-m",   
        dest = "mi_in", 
        help = "[Required] mi summary file")
    parser.add_option("-o", 
        dest = "mi_out", 
        help = "[Required] filtered mi summary file")


    filtering_options = OptionGroup(parser, "Filtering Options")

    filtering_options.add_option("--AR_min", 
        dest = "AR_min",     
        default = 0.1, type = "float", 
        help = "Minimum allelic ratio of Haplotype [Default: 0.1]")
    filtering_options.add_option("--AR_max", 
        dest = "AR_max",     
        default = 0.9, type = "float", 
        help = "Maximum allelic ratio of Haplotype [Default: 0.9]")
    filtering_options.add_option("--PSI_min", 
        dest = "PSI_min",     
        default = 0.1, type = "float", 
        help = "Minimum PSI of exonic part [Default: 0.1]")
    filtering_options.add_option("--PSI_max", 
        dest = "PSI_max",     
        default = 0.9, type = "float", 
        help = "Maximum PSI of exonic part [Default: 0.9]")
    filtering_options.add_option("--COV_min", 
        dest = "COV_min",     
        default = 20, type = "int", 
        help = "Minimum total coverage per phasing group [Default: 20]")
    filtering_options.add_option("--delta_PSI_min", 
        dest = "delta_PSI_min",     
        default = 0.01, type = "float", 
        help = "Minimum delta PSI between haplotypes [Default: 0.05]")
    filtering_options.add_option("--exonic_part_len_min", 
        dest = "exonic_part_len_min",     
        default = 3, type = "int", 
        help = "Minimum length of exonic part [Default: 3]")
    filtering_options.add_option("--linkage_types", 
        dest = "linkage_types",     
        default = "ASTS,HAP-ExonicPart", type = "str", 
        help = "Coma-separated list of labels to keep [Default: 'ASTS,HAP-ExonicPart']")

    parser.add_option_group(filtering_options)

    
    (options, args) = parser.parse_args()

    acceptable_linkage_types  = options.linkage_types.split(',') 

    df = pd.read_csv(options.mi_in, sep = '\t', escapechar = '#')

    df = df[(df['label'] == "cis")]
    logging.info(f"Number of cis events: {len(df)}")

    df = df[(df['Linkage_Type'].isin( acceptable_linkage_types ))]
    logging.info(f"Number of right linakage types: {len(df)}")

    df = df[(df['Cov'] >= options.COV_min)]
    logging.info(f"Number of events with good coverage: {len(df)}")

    df = df[((df['end'] - df['start']) >= options.exonic_part_len_min)]
    logging.info(f"Number of events with long exonic parts: {len(df)}")

    df = df[(df['PSI'] >= options.PSI_min)]
    df = df[(df['PSI'] <= options.PSI_max)]
    logging.info(f"Number of events within PSI ratio: {len(df)}")

    df = df[(df['AR']  >= options.AR_min)]
    df = df[(df['AR']  <= options.AR_max)]
    logging.info(f"Number of events within AR ratio: {len(df)}")

    df = df[(abs(df['Delta_PSI']) >= options.delta_PSI_min)]
    logging.info(f"Number of events with high delta PSI: {len(df)}")
   
    df.to_csv(options.mi_out, sep = '\t', index = False)


if __name__ == "__main__":
    main()
