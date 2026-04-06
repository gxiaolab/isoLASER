#!/usr/bin/python
import sys
import re
import sklearn.metrics
import HTSeq
import pysam
import pickle
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from scipy.stats import chisquare, norm, chi2_contingency

from . import gqv_software_management as SOFTMAN


def split_GenomicInterval_string(gi_string):
    parsed_groups = re.search(r'(\w+):\[(\d+),(\d+)\)/([+-])', gi_string)
    chrom, start, end, strand = parsed_groups.groups()
    return (chrom, int(start), int(end), strand)



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



def get_genes_from_gtf(gff_file, region_chrom = None, region_start = None, region_end = None):

    gene_coords_dict = dict() 

    df = pd.read_csv(f'{gff_file}/GeneList.tsv', sep = '\t', header = None)
    df.columns = ['gene_id', 'gene_name', 'gi_interval']

    for row in df.itertuples():

        chrom, start, end, strand = split_GenomicInterval_string(row.gi_interval)

        gene_name = f"{row.gene_id}_{row.gene_name}"

        if region_chrom is None:
            gene_coords_dict[gene_name] = (chrom, start, end, strand)
        elif chrom == region_chrom and region_start is None and region_end is None:
            gene_coords_dict[gene_name] = (chrom, start, end, strand)      
        elif chrom == region_chrom and start < region_end and end > region_start:
            gene_coords_dict[gene_name] = (chrom, start, end, strand)

    return gene_coords_dict 



def mutual_information(hap_clusters, intron_reads,
    min_comb_cov = 6, part_alleles = []):
    
    ## Making matrix 
    common_reads = list(set(hap_clusters) & set(intron_reads))
    cov = len(common_reads)

    mat = np.zeros((cov, 2)).astype(int)

    for read_i, r in enumerate(common_reads):
        mat[read_i][0] = hap_clusters[r] + 1 # (1, 2) 
        mat[read_i][1] = intron_reads[r] + 1 # (1, 2)

    unique_ip = np.unique(mat[:, 1])
   
    ## Mat format 
    mat_u = np.zeros((2, 2)).astype(int)
    mat_combs, mat_counts = np.unique(mat, axis = 0, return_counts = True)
    for (_hap, _ip), count in zip(mat_combs, mat_counts):
        mat_u[_hap - 1][_ip - 1] = count 

    ## Filter based on coverage 

    hap_covs = np.sum(mat_u, axis = 1) # rowsums

    A_freq = mat_u[0]/np.sum(mat_u[0])
    B_freq = mat_u[1]/np.sum(mat_u[1])

    if np.all(hap_covs >= min_comb_cov): #and np.all(colSums >= min_comb_cov):
        AMI_score = sklearn.metrics.adjusted_mutual_info_score(mat[:, 0], mat[:, 1]) * np.log(2)
        delta_psi = abs(A_freq[1] - B_freq[1])
    else:
        AMI_score = np.nan
        delta_psi = np.nan

    ## String format 

    linkage_string = []

    for _hap in [1, 2]:
        for _ip in [1, 2]:
            linkage_string.append( mat_u[_hap - 1][_ip - 1] )
    
    linkage_string = ';'.join(map(str, linkage_string))

    return AMI_score, cov, linkage_string, delta_psi



def process_reads_into_txs(RCG_reads, TX_structure):

    RCG_transcripts = defaultdict(set)
    for read_name, attr in RCG_reads.items():
        transcript_id = attr['ZT']
        RCG_transcripts[transcript_id].add(read_name)

    IntronReads = defaultdict(dict) 
    

    for (part_annot, ExonicPart), part_alleles in TX_structure.items():

        if part_annot != "exonic_part":
            continue 

        part_alleles_sorted = sorted(list(part_alleles.keys()))
        
        for j, exon_or_intron in [(0, 'intron'), (1, 'exon')]: 
            for transcript_id in part_alleles[exon_or_intron]:
                for read_name in RCG_transcripts[transcript_id]:
                    IntronReads[(part_annot, ExonicPart)][read_name] = j

    return IntronReads



def asts(RCG_reads, hap_clusters):

    MAT = defaultdict(lambda : defaultdict(lambda : 0))

    for read_name, attr in RCG_reads.items():

        transcript_id = attr['ZT']
        if read_name not in hap_clusters:
            continue 
        cluster_id = hap_clusters[read_name]
        MAT[cluster_id][transcript_id] += 1

    MAT = pd.DataFrame(MAT)
    MAT = MAT.fillna(0)

    if MAT.shape[1] == 1:
        return np.nan, np.nan, np.nan
    else:
        row_sums = MAT.sum(axis = 1)
        MAT = MAT[row_sums > 0]

        real_counts = [list(MAT[0]), list(MAT[1])]
        chi2, p_value, a, b = chi2_contingency(real_counts)

        allelic_string = []
        MAT = MAT.loc[row_sums.sort_values(ascending = False).index]

        for row in MAT.itertuples():
            allelic_string.append(f'{row._1:.0f},{row._2:.0f}')
        allelic_string = ":".join(allelic_string)

        return p_value, 'na', allelic_string



def mi_parse_variants(ALL_VARS_FEAT, RCG_reads, gene_pickle_file, RC_info, all_clusters, polyfits, normalfits):


    (Gene_name, CHROM, GENE_Start, GENE_End, STRAND), RCB_list = RC_info
    Gene_coords = f'{CHROM}\t{GENE_Start}\t{GENE_End}\t{STRAND}'

    with open(gene_pickle_file, 'rb') as tb:
        TX_structure = pickle.load(tb)

    IntronReads = process_reads_into_txs(RCG_reads, TX_structure)

    printing_lines = []
    haplotagged_reads = defaultdict(dict)

    ### Haplotype - exonic_part linkage

    for phasing_group, clusters in all_clusters:
        
        hap_clusters = {}

        for cluster_id, read_dict in clusters.items():
            for read_i in read_dict:
                hap_clusters[read_i] = cluster_id 
                haplotagged_reads[(read_i, phasing_group)] = (cluster_id, Gene_name, CHROM)

        ## ASTS 
        asts_chi_p, label, asts_allelic_string = asts(RCG_reads, hap_clusters)
        
        hap_count = Counter(hap_clusters.values())
        coverage  = sum(hap_count.values())

        allelic_ratio = [v/coverage for v in hap_count.values()][0]
        allelic_ratio = round(allelic_ratio, 3)

        outline1 = f"{Gene_coords}\t{Gene_name}\tASTS\tHAP-{phasing_group}"
        outline2 = f"{np.nan}\t{asts_chi_p:.2e}\t{coverage}\t{np.nan}\t{allelic_ratio:.3f}\t{asts_allelic_string}\t{label}\t{np.nan}"

        printing_lines.append(outline1 + '\t' + outline2)

        ### Label RNA editing sites gf

        for POS, ALLELE_dict in ALL_VARS_FEAT.items():
            for ALLELE, attr in ALLELE_dict.items():

                if not attr["FEAT"]['is_editing']:
                    continue

                rnae_reads = attr["READS"]            
                
                mi_out = mutual_information(hap_clusters, rnae_reads)

                (ami, lc, link_str, delta_psi) = mi_out

                editing_ratio = np.mean(list(rnae_reads.values()))

                if ami < 0.01 and editing_ratio > 0.05 and editing_ratio < 0.95 and sum(rnae_reads.values()) >= 3:
                    ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["variant_type"] = "EditingSite"
                else:
                    ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["variant_type"] = "SNV"

        ### Haplotype - exonic_part linkage
        
        for (part_annot, ExonicPart), intron_reads in IntronReads.items():
            
            psi = np.array(list(intron_reads.values()))
            psi = np.mean(psi) 
            
            mi_out = mutual_information(hap_clusters, intron_reads)

            (ami, cov, allelic_string, delta_psi) = mi_out 

            label, norm_p = get_ami_cutoffs(polyfits, normalfits, cov, ami)  

            ExonicCoord = '\t'.join(map(str, split_GenomicInterval_string(ExonicPart)))       
           
            outline1 = f"{ExonicCoord}\t{Gene_name}\tHAP-ExonicPart\tHAP-{phasing_group}"
            outline2 = f"{ami:.3f}\t{norm_p:.3e}\t{cov}\t{psi:.3f}\t{allelic_ratio:.3f}\t{allelic_string}\t{label}\t{delta_psi:.3f}"

            printing_lines.append(outline1 + '\t' + outline2)

           
    for POS, ALLELE_dict in ALL_VARS_FEAT.items():
        for ALLELE, attr in ALLELE_dict.items():

            var_reads = attr["READS"]          
            var_id    = attr['FEAT']['variant_id']
            var_gt    = attr['FEAT']['GT']
            
            allelic_ratio = np.nanmean(list(var_reads.values()))

            for (part_annot, ExonicPart), intron_reads in IntronReads.items():
            
                psi = np.array(list(intron_reads.values()))
                psi = np.mean(psi) 

                mi_out = mutual_information(var_reads, intron_reads)

                (ami, cov, allelic_string, delta_psi) = mi_out 

                if cov < 10:
                    continue

                ExonicCoord = '\t'.join(map(str, split_GenomicInterval_string(ExonicPart))) 

                label, norm_p = get_ami_cutoffs(polyfits, normalfits, cov, ami)                            
                outline1 = f"{ExonicCoord}\t{Gene_name}\tVAR\t{var_id}"
                outline2 = f"{ami:.3f}\t{norm_p:.3e}\t{cov}\t{psi:.3f}\t{allelic_ratio:.3f}\t{allelic_string}\t{label}\t{delta_psi:.3f}"
                printing_lines.append(outline1 + '\t' + outline2)

    return printing_lines , haplotagged_reads
