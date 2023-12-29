#!/usr/bin/python
import sys
import subprocess
import pandas as pd
from optparse import OptionParser, OptionGroup
import os


def test_gatk():
    try:
        subprocess.run(['gatk', '--java-options', '-Xmx4g', '--version'], check = True)
    except:
        print("ERROR: GATK not found. Please install GATK and add it to your PATH.")
        sys.exit(1)

def main():

    usage = f"\n\tpython {sys.argv[0]} -f <fofn.txt> -o <outfile.txt>"

    parser = OptionParser(usage = usage, description = '')

    parser.add_option("-i", dest = "fofn")
    parser.add_option("-o", dest = "prefix")
    parser.add_option("-f", dest = "REF_file")

    (options, args) = parser.parse_args()

    test_gatk()

    fofn_file = f"{options.prefix}.gvcfs.list"
    vcf_file  = f"{options.prefix}.merged.gvcf.gz"
    gvcf_file = f"{options.prefix}.merged.genotyped.gvcf.gz"

    df = pd.read_table(options.fofn, header = None, names = ["sm", "mi", "vcf"], sep = "\t") 
    df[["vcf"]].to_csv(fofn_file, header = False, index = False, sep = " ")

    ## Combine gvcfs

    CMD_LINE_GATK_CombineGVCFs = ['gatk', '--java-options', '-Xmx4g', 'CombineGVCFs',
                                    '--variant', fofn_file, 
                                    '-R', options.REF_file, 
                                    '-O', vcf_file]

    subprocess.run(CMD_LINE_GATK_CombineGVCFs, check = True)

    ## Genotype gvcfs

    CMD_LINE_GATK_GenotypeGVCFs = ['gatk', '--java-options', '-Xmx4g', 'GenotypeGVCFs',
                                    '-R', options.REF_file, 
                                    '-V', vcf_file, 
                                    '-O', gvcf_file]

    subprocess.run(CMD_LINE_GATK_GenotypeGVCFs, check = True)

    ## Index file

    CMD_LINE_TABIX_Index = ['tabix', '-p', 'vcf', '-f', gvcf_file]

    subprocess.run(CMD_LINE_TABIX_Index, check = True)


if __name__ == "__main__":
    main()