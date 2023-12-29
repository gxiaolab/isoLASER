#!/usr/bin/python
import numpy as np
from optparse import OptionParser
from collections import defaultdict
from scipy.stats import chisquare, norm

import vcf
import sys
import os
import re
import pickle
import sklearn.metrics
import pandas as pd 
import tabix
import time
import HTSeq
import copy
import warnings
import logging

from . import gqv_gff_utils 
warnings.filterwarnings("ignore")



def print_time_stamp(message = ""):
    sys.stderr.write(message + '\n')
    sys.stderr.flush()


def get_sqtl_slope(data):
    data = np.array(data)
    try:
        slope, offset = np.polyfit(data[:, 0], data[:, 1], 1, w = np.log10(data[:,2] + 1))
        return round(slope, 3)
    except:
        return np.nan


def get_ami_cutoffs(polyfits, normalfits, cov, test_ami):
    log10cov = np.log10(cov)
    label = 'trans' 

    for quantile, ami_label in [(0.95, 'amb'), (0.99, 'cis')]:
        q_dict = polyfits[quantile]
        amis = [func(log10cov) for func in q_dict.values()]
        ami_mean = np.mean(amis)
        
        if test_ami > ami_mean: 
            label = ami_label 
        else:
            break 

    if np.isnan(test_ami):
        label = np.nan

    min_cov = normalfits[0][0]
    max_cov = normalfits[-1][0]
    N = len(normalfits)
    log10cov = round(log10cov, 3)

    if log10cov <= min_cov:
        idx = 0
        log10cov = min_cov  
    elif log10cov >= max_cov:
        idx = N - 1 
        log10cov = max_cov
    else:
        idx = int((log10cov - min_cov) * N/(max_cov - min_cov))

    p = np.nan
    for row in normalfits[max(idx - 3 , 0) : idx + 3]:
        mu_cov, std_cov, mu_ami, std_ami, mu_delta_psi, std_delta_psi = row
        mu_cov, std_cov = round(mu_cov, 3), round(std_cov, 3)

        if log10cov >= mu_cov - std_cov and log10cov <= mu_cov + std_cov:
            p = 1 - norm.cdf(test_ami, mu_ami, std_ami)
            break

    return label, p 


def calculate_ami(counts_dict, min_n_obs = 10):
    
    c = [[counts_dict["REF_e"], counts_dict["REF_i"]],\
         [counts_dict["ALT_e"], counts_dict["ALT_i"]]]

    c = np.array(c)
    HapCov  = np.sum(c, axis = 1)
    ExonCov = np.sum(c, axis = 0)

    if not np.all(HapCov >= min_n_obs):
        return np.nan, np.nan, np.nan, np.nan

    mat = []
    for k, count in counts_dict.items():
        allele, part = k.split("_")
        mat += ([[allele, part]] * count)

    mat = np.array(mat)
    AMI_score = sklearn.metrics.adjusted_mutual_info_score(mat[:, 0], mat[:, 1]) * np.log(2)

    psi = c[:, 0]/HapCov  
    delta_psi = psi[0] - psi[1]
    combined_psi = ExonCov[0]/np.sum(ExonCov)
    combined_ar  = HapCov[0]/np.sum(HapCov)

    return round(AMI_score, 3), round(delta_psi, 3), round(combined_psi, 3), round(combined_ar, 3)



def read_mi(fofn_df, Gene_Id, Gene_attr):    

    exonic_parts = defaultdict(lambda : defaultdict(dict))
    CHROM, START, END, STRAND = Gene_attr

    for metadata_row in fofn_df.itertuples():


        mi_tabix = tabix.open(metadata_row.mi)
        header = pd.read_csv(metadata_row.mi, nrows = 0, sep = '\t').columns
        records = list(mi_tabix.query(CHROM, int(START), int(END)))

        mi_subset = pd.DataFrame(records, columns = header) 
        mi_subset = mi_subset[(mi_subset['Gene_name'] == Gene_Id)]
        mi_subset = mi_subset[(mi_subset['Linkage_Type'] == 'VAR')]

        for row in  mi_subset.itertuples():
        
            d_keys = ['REF_e', 'REF_i', 'ALT_e', 'ALT_i']
            d_vals = map(int, row.allelic_counts.split(';'))
            hap_counts = dict(zip(d_keys, d_vals))

            ExPart = f"{CHROM}:{row.start}-{row.end}:{row.strand}"

            DP  = sum(hap_counts.values())
            AC  = hap_counts["ALT_e"] + hap_counts["ALT_i"]
            AB  = round((AC + 0.05)/(DP + 0.05), 3)
            
            PSI = (hap_counts["REF_i"] + hap_counts["ALT_i"] + 0.05)/(DP + 0.05)
            PSI = round(PSI, 3)

            exonic_parts[row.PG][metadata_row.sm][ExPart] = (AB, PSI, float(row.AMI), DP, hap_counts)
    
    return exonic_parts 



def process_variant(exonic_parts, VCF_handle):
 
    for variant in list(exonic_parts):
        if ">" in variant:
            chrom, pos, alleles = variant.split(':')
            found = False
            try:
                for record in VCF_handle.fetch(chrom, int(pos) - 3, int(pos) + 3):
                    if record.POS - 1 == int(pos):
                        found = True
            except:
                pass 

            if not found:
                exonic_parts.pop(variant)


    for ref_block in list(exonic_parts):
        if ">" not in ref_block: 
            tmp_dict = exonic_parts.pop(ref_block)

            r = re.match('(chr\w+):\[(\d+),(\d+)\)(.)', ref_block)

            chrom = r.group(1)
            REF_block_start = int(r.group(2)) # 0-base
            REF_block_end   = int(r.group(3)) # 1-base 
            
            try: 
                for record in VCF_handle.fetch(chrom, REF_block_start - 1, REF_block_end + 1):
                    real_variant = f"{chrom}:{record.POS - 1}:{record.REF}>{record.ALT[0]}" 

                    if real_variant not in exonic_parts:
                        continue
                    
                    record_POS = record.POS - 1 #0-base
                    record_END = record.POS - 1 + len(record.REF) #1-base

                    if REF_block_start <= record_POS and REF_block_end >= record_END:
                        
                        for sample_name in tmp_dict:

                            var_sample_info = [s for s in record.samples if s.sample == sample_name][0]
                            
                            for part in tmp_dict[sample_name]:
                                vals = copy.deepcopy( tmp_dict[sample_name][part] )

                                if part in exonic_parts[real_variant][sample_name]:

                                    ab_diff = abs(var_sample_info['AB'][0] - vals[0])

                                    if ab_diff < 0.05 and vals[3] > exonic_parts[real_variant][sample_name][part][3]:

                                        exonic_parts[real_variant][sample_name][part] = vals

                                elif part not in exonic_parts[real_variant][sample_name]:
                                    exonic_parts[real_variant][sample_name][part] = vals
                  
            except Exception as e:
                sys.stderr.write(f'Exception {e}\n')
                pass 


def norm_cov(HAP_COUNTS):

    covs = [sum(h.values()) for h in HAP_COUNTS]
    
    min_cov = np.mean(covs)

    HAP_COUNTS_FINAL = defaultdict(lambda : 0)

    for hap_counts in HAP_COUNTS:
        cov = sum(hap_counts.values())
        scaling_ratio = min_cov/cov  
        for hap, N in hap_counts.items():
            new_count = int(round(scaling_ratio * N, 0))
            HAP_COUNTS_FINAL[hap] += new_count

    return HAP_COUNTS_FINAL



def merge_mi_summary_files(fofn_df, vcf_file, out_file, polyfits, normalfits, gene_lst):

    SMList = list(fofn_df.sm.unique())

    ## output file
    outf = open(out_file, 'w')

    # #chrom	start	end	strand	Gene_name	Linkage_Type	PG	AMI	chi_p	Cov	PSI	AR	H1e;H1i;H2e;H2i	label	Delta_PSI
    out_line = "Gene_name\tvariant_qual\texonic_part\tcombined_AMI\tlabel\tdelta_psi(ref-alt)\tCOV\tFORMAT\t" 
    out_line += "\t".join(SMList) + "\tSlope\tPSI\tAR\tp_val"
    outf.write(out_line + '\n')

    VCF_handle = vcf.Reader(filename = vcf_file)

    for g_i, (Gene, Gene_coords) in enumerate(gene_lst.items()):

        if g_i % 1000 == 0: 
            logging.info(f"{g_i:,} genes processed\n")

        ExonicPartsDict = read_mi(fofn_df, Gene, Gene_coords)

        process_variant(ExonicPartsDict, VCF_handle)
        
        for Variant, BamDict in ExonicPartsDict.items():

            regions = defaultdict(set)
            for Bam in SMList:
                for region in BamDict[Bam]:
                    regions[region].add(Bam)

            for ExPart in regions:
                
                if len(regions[ExPart]) < 2:
                    continue 

                values = [] 
                HAP_COUNTS_tmp = []
                slope_data = []

                for Bam in SMList: # (AB, PSI, AMI, AMI, hap_counts)
                    try:
                        vals = BamDict[Bam][ExPart]
                        AB, PSI, AMI, DP, hap_counts = vals 
                    except KeyError:
                        AB, PSI, AMI, DP, hap_counts = np.nan, np.nan, np.nan, 0, {}
                    
                    values.append(f"{AB:.2f}:{PSI:.2f}:{AMI:.2f}:{DP}")
                    slope_data.append([AB, PSI, DP])
                    
                    if DP >= 10:
                        HAP_COUNTS_tmp.append(hap_counts)

                # norm read count
                
                HAP_COUNTS = norm_cov(HAP_COUNTS_tmp)

                combined_cov = sum(HAP_COUNTS.values()) 
                combined_ami, delta_psi, combined_psi, combined_ar = calculate_ami(HAP_COUNTS)
                slope = get_sqtl_slope(slope_data)

                if np.isnan(combined_ami):
                    # print('skipping', Gene, Variant, ExPart, HAP_COUNTS)
                    continue 

                label, p = get_ami_cutoffs(polyfits, normalfits, combined_cov, combined_ami)
                vals_string = "AB:PSI:AMI:DP\t" + "\t".join(values)

                outline = [Gene, Variant, ExPart, combined_ami, label, delta_psi, combined_cov, vals_string, slope, combined_psi, combined_ar, p]
                outf.write('\t'.join(map(str, outline)) + '\n')
        
    outf.close()
        




def main():
    usage = f"\n\tpython {sys.argv[0]} -f <fofn.txt> -o <outfile.txt>"

    sys.stderr.write("CMD = python" + " ".join(sys.argv) + "\n")

    parser = OptionParser(usage = usage, description = '')

    parser.add_option("-i", dest = "fofn")
    parser.add_option("-o", dest = "outf")
    parser.add_option("-v", dest = "vcf")
    parser.add_option("-t", dest = "transcript_db")

    (options, args) = parser.parse_args()

    PACKAGE_PATH = os.path.realpath(os.path.dirname(__file__))

    with open(f'{PACKAGE_PATH}/glm/ami_null.polyfits.ar_psi_bins.pickle', 'rb') as f:
        polyfits = pickle.load(f)
    with open(f'{PACKAGE_PATH}/glm/mi_testing.training_summary.normalfits.pickle', 'rb') as f:
        normalfits = pickle.load(f)

        
    gene_lst = gqv_gff_utils.get_genes_from_gtf(options.transcript_db)

    df = pd.read_table(options.fofn, header = None, names = ["sm", "mi", "vcf"], sep = "\t") 

    merge_mi_summary_files(df, options.vcf, options.outf, polyfits, normalfits, gene_lst)
    sys.stderr.write(f"\t[Progress] done\n" ); sys.stderr.flush()


if __name__ == "__main__":
    main()