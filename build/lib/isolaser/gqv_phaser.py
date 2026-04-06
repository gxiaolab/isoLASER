#!/usr/bin/python
import sys
import numpy as np
import scipy.sparse
import copy
import random
import itertools
from collections import defaultdict, Counter
from sklearn.cluster import SpectralClustering
from time import time
import warnings

from . import gqv_software_management as SOFTMAN
from . import gqv_glm


GENOTYPE_OPTIONS = [[0, 0], [0, 1], [1, 1]]
GENOTYPE_OPTIONS_EXT = [[0, 0], [0, 1], [1, 1], [1, 2], [2, 2]]




def main(OUTPUT_LINES, RCG_reads, clf_dict, options):
    """
    time every step of this function
    """
    RCG_variants, RCG_reads, mat = table_init(OUTPUT_LINES, RCG_reads)
    
    for _ in range(3):
        QUAL_weigths = update_pl_gt(OUTPUT_LINES, RCG_variants, mat)

        all_clusters = table_genotyping(OUTPUT_LINES, RCG_variants, RCG_reads, mat, QUAL_weigths, options, _)

        gqv_glm.update_vars(clf_dict, OUTPUT_LINES, options)

    refine_variants(OUTPUT_LINES, RCG_variants, mat)   

    return OUTPUT_LINES, all_clusters



def table_init(OUTPUT_LINES, RCG_reads):

    RCG_variants = []

    for anchor_pos, anchor_alleles in OUTPUT_LINES.items():
        for allele in anchor_alleles:
            RCG_variants.append([anchor_pos, allele])

    RCG_variants.sort(key = lambda h : h[0][2])
    RCG_reads = np.array(RCG_reads)

    n = len(RCG_variants)
    m = len(RCG_reads)

    mat = np.zeros((m, n))
    
    if n >= 1 and m > 1:     
        for j, (POS, ALLELE) in enumerate(RCG_variants):
            for i, read_name in enumerate(RCG_reads):
                var_reads = OUTPUT_LINES[POS][ALLELE]["READS"]
                if read_name in var_reads:
                    mat[i][j] = 1 + var_reads[read_name] 

    mat[mat == 0] = np.nan
    
    return RCG_variants, RCG_reads, mat



def update_pl_gt(OUTPUT_LINES, RCG_variants, mat):
    QUAL_weigths = []

    for j, (POS, ALLELE) in enumerate(RCG_variants):
        
        feat = OUTPUT_LINES[POS][ALLELE]['FEAT']

        QUAL = feat["QUAL"]
        QUAL_weigths.append(QUAL)

        PL = PL_prob(mat[:, j] - 1, [1.0 , max(1e-3, QUAL)] )
        feat['PL'] = PL
                                                            
        GT_PL = GENOTYPE_OPTIONS[np.argmin(PL)]
        feat['GT'] = GT_PL 

        try:
            PL_set = sorted(set(PL))
            GQ = min(PL_set[1] - PL_set[0], 99)
        except:
            GQ = 1
        
        feat['GQ'] = GQ

    return QUAL_weigths




def table_genotyping(OUTPUT_LINES, RCG_variants, RCG_reads, mat, QUAL_weights, options, turn):
    all_clusters = []

    ### Use heterozygous SNPs for phasing 0.05 <= AR <= 0.95  
    vars_ab = np.nanmean(mat, axis = 0)
    vars_type = np.array([rc[0][3] for rc in RCG_variants])
    vars_snps = vars_type == "SNV"

    vars_hete = (vars_ab >= 1.05) & (vars_ab <= 1.95) 
    vars_homo = (vars_ab < 1.05)  | (vars_ab > 1.95)

    snps_hete = vars_hete * vars_snps 

    PhasGrp = 0

    for var_type, ploidy in [(snps_hete, 2), (vars_homo, 1)]:

        dmat = np.copy(mat)
    
        connected_reads_subgraph, connected_variants_subgraph = find_connected_components(dmat, var_type, min_rn = options.minCoverage)
        
        for sub_ReadList, sub_VarList in zip(connected_reads_subgraph, connected_variants_subgraph):

            if len(sub_ReadList) < options.minCoverage:
                continue

            sub_mat = mat[list(sub_ReadList), ] 
            sub_rns = RCG_reads[list(sub_ReadList)]
           
            weights = np.array(QUAL_weights)
            haplotypes, clusters = k_means_clustering(sub_mat, weights, ploidy = ploidy)

            hap_mat = np.zeros(sub_mat.shape) 
            hap_mat[:, :] = 1.5 

            for hap_i, read_dict in list(clusters.items()):
                hap_mat[list(read_dict)] = haplotypes[hap_i] 
                clusters[hap_i] = [sub_rns[ri] for ri in read_dict]

            all_clusters.append((PhasGrp, clusters))

            diff_mat = np.abs(hap_mat - sub_mat) 

            for j, (POS, ALLELE) in enumerate(RCG_variants):
                if j in sub_VarList:
                    if (vars_hete[j] and ploidy == 2) or (vars_homo[j]):

                        normHapDistance = np.log(np.nanmean(diff_mat[:, j]) + 1e-5)
                        OUTPUT_LINES[POS][ALLELE]['FEAT']['ndHAP'] = normHapDistance    
                                                                            
                        GT_array = haplotypes[:, j] - 1 
                        if ploidy == 1:
                            GT_array = [GT_array[0], GT_array[0]]

                        OUTPUT_LINES[POS][ALLELE]['FEAT']['GT'] = GT_array 
                        OUTPUT_LINES[POS][ALLELE]['FEAT']['PG'] = PhasGrp

            PhasGrp += 1

    return all_clusters
    


def PL_prob(r_array, Q, genotype_options = GENOTYPE_OPTIONS):
    PL_array = []
    r_array = np.array(r_array)

    for (gt1, gt2) in genotype_options:
        Q_gt1 = Q[gt1]
        Q_gt2 = Q[gt2]
        
        p = np.full_like(r_array, 1e-5)
        p += 0.5 * ((r_array == gt1).astype(int) * Q_gt1 + (1 - (r_array == gt1).astype(int)) * (1 - Q_gt1))
        p += 0.5 * ((r_array == gt2).astype(int) * Q_gt2 + (1 - (r_array == gt2).astype(int)) * (1 - Q_gt2))
        
        PL = np.sum(np.log10(p))
        PL_array.append(PL)
        
     
    PL_array = np.array(PL_array) * -10 
    PL_array = PL_array - np.min(PL_array)
    
    return PL_array.astype(int)





# def PL_prob(r_array, Q, genotype_options = GENOTYPE_OPTIONS):
#     PL_array = []
#     for (gt1, gt2) in genotype_options:
#         PL = 0.0 
#         for r in r_array:

#             p = 1e-5
#             p += 0.5 * ((r == gt1)*Q[gt1] + (1 - (r == gt1))*(1 - Q[gt1]))
#             p += 0.5 * ((r == gt2)*Q[gt2] + (1 - (r == gt2))*(1 - Q[gt2]))

#             PL += np.log10(p)       
        
#         PL_array.append(PL)

#     PL_array = [int(-10 * pl) for pl in PL_array]
#     PL_array = [(pl - min(PL_array)) for pl in PL_array]

#     return PL_array





def find_connected_components(dmat, subsetting_array, min_rn = 10):
    '''
    dmat = sparse matrix (reads x variants)
    dmat1 = sparse matrix (variants x variants)
    ''' 
    dmat[dmat >= 1] = 1
    dmat[np.isnan(dmat)] = 0

    dmat = dmat * subsetting_array
    dmat = scipy.sparse.coo_matrix(dmat, dtype = int)
    dmat1 = dmat.T.dot(dmat) 

    ## variants with at least min_rn reads can be phased together
    dmat1[dmat1 < min_rn] = 0
    dmat1.eliminate_zeros() 

    visited = {}
    for i in range(len(subsetting_array)):
        if subsetting_array[i]:
            visited[i] = False 

    connected_variants_subgraph = []

    def dfs(mat_i):
        for mat_j in dmat1[mat_i].indices:
            if mat_i == mat_j:
                continue
            if not visited[mat_j]:
                visited[mat_j] = True 
                dfs(mat_j)

    prev_set = set()
    for i in range(dmat1.shape[0]):
        if i not in visited:
            continue
        if not visited[i]:
            visited[i] = True 
            dfs(i)
            connected_set = set([k for k, v in visited.items() if v]) - prev_set
            connected_variants_subgraph.append(connected_set)
            prev_set |= connected_set


    connected_reads_subgraph = []
    for connected_nodes in connected_variants_subgraph:
        read_set = set()
        for vi in connected_nodes:
            read_set |= set(dmat.row[dmat.col == vi])
        connected_reads_subgraph.append(read_set)

    return connected_reads_subgraph, connected_variants_subgraph


def vector_dist(v1, v2, w):
    raw  = (v1 - v2)*w 
    dist = np.sqrt(np.nansum(np.square(raw)))
    return dist



def k_means_clustering(mat, w, ploidy, iterations = 5, random_initialization = 10):

    m, n = mat.shape   
    score_track = {}

    non_zero_weights = np.nansum(w > 0)

    for _rinit in range(min(random_initialization, non_zero_weights + 1)):

        centroids = mat[np.random.randint(0, m, ploidy), :]
        prev_centroids = centroids.copy()

        for _iter in range(iterations):
            clusters = defaultdict(dict)
            
            for i in range(m):
                cl_dist = {}
                for PhasGrp in range(ploidy):
                    cl_dist[PhasGrp] = vector_dist(centroids[PhasGrp], mat[i], w) 

                best_centroid = min(cl_dist, key = cl_dist.get)
                clusters[best_centroid][i] = cl_dist[best_centroid]

            read_dist2hap = []

            for PhasGrp, rows in clusters.items():
                tmp = np.nanmean(mat[list(rows), :], axis = 0)
                tmp[np.isnan(tmp)] = centroids[PhasGrp][np.isnan(tmp)]

                centroids[PhasGrp] = tmp

                read_dist2hap.extend(rows.values())

            if (centroids == prev_centroids).all():
                break
            else:
                prev_centroids = centroids.copy()

        for PhasGrp, rows in clusters.items():
            centroids[PhasGrp] = np.nanmean(mat[list(rows), :], axis = 0)

        score_track[(_rinit, ploidy)] = [np.sum(read_dist2hap), np.round(centroids, 0), clusters]

    best_r_i = min(score_track, key = lambda h : score_track[h][0])
    best_v, centroids, clusters = score_track[best_r_i]
    
    return centroids, clusters


def refine_variants(OUTPUT_LINES, RCG_variants, mat):

    allele_groups = defaultdict(list)

    for j, (genomic_pos, genomic_allele) in enumerate(RCG_variants):
        feat = OUTPUT_LINES[genomic_pos][genomic_allele]['FEAT']
        qual = feat['QUAL']
        gt   = feat['GT']
        pg   = feat['PG']

        (CHROM, Gene_Id, POS, var_group) = genomic_pos

        if var_group == "REF_BLOCK":
            feat['GT_str'] = "0/0"
            continue

        allele_groups[genomic_pos].append((genomic_allele, j, qual, gt, pg))
        # print('\ttmp', genomic_pos, genomic_allele, j, qual, gt, pg)

    for genomic_pos, genomic_allele_list in allele_groups.items():
        sub_mat = np.copy(mat) 

        genomic_allele_list.sort(key = lambda x : x[2], reverse = True)
        
        for val in genomic_allele_list[:]:
            allele, j, qual, gt, pg = val 
            if np.sum(gt) == 0.0 or qual < 0.5:
                genomic_allele_list.remove(val)

        for i, val in enumerate(genomic_allele_list[:]):
            allele, j, qual, gt, pg = val 
            sub_mat[:, j][sub_mat[:, j] == 2] += i 
            if i >= 2: 
                genomic_allele_list.remove(val)


        sub_mat = sub_mat[: , [v[1] for v  in genomic_allele_list]]
        sub_mat = sub_mat - 1.0

        row_covs = np.nansum(sub_mat > 0, axis = 1)      
        if not np.all(row_covs <= 1):
            print(  f"error {genomic_pos} \t {list(row_covs)},{ [v[1] for v  in genomic_allele_list] }"  )         
        
        idx_na = np.all(np.isnan(sub_mat), axis=1)

        row_sums = np.nansum(sub_mat, axis = 1)
        row_sums = row_sums[~idx_na]

        # genotyeps pl
        Q = [1.0] + [v[2] for v  in genomic_allele_list]
        Q = [max(1e-3, q) for q in Q]
        
        if len(Q) == 1:
            pl = PL_prob(row_sums, Q + [np.nan])
        elif len(Q) == 2:
            pl = PL_prob(row_sums, Q)
        else:
            pl = PL_prob(row_sums, Q, genotype_options = GENOTYPE_OPTIONS_EXT)
        
        # most likely genotype
        GT_PL = GENOTYPE_OPTIONS_EXT[np.argmin(pl)]


        # gq
        try:
            PL_set = sorted(set(pl))
            GQ = min(PL_set[1] - PL_set[0], 99)
        except:
            GQ = 1

        # phasing
        # if all(v[4] != 'NA' for v  in genomic_allele_list) and GT_PL[0] != GT_PL[1]:
        #     GT_str = "{0}|{1}".format(*GT_PL)
        # else:
        #     GT_str = "{0}/{1}".format(*GT_PL)


        for val in genomic_allele_list:
            genomic_allele, j, qual, gt, pg = val 

            if pg != "NA" and gt[0] != gt[1] and sum(gt) == sum(GT_PL):
                GT_str = "{0:.0f}|{1:.0f}".format(*gt)
            else:
                GT_str = "{0}/{1}".format(*GT_PL)

            OUTPUT_LINES[genomic_pos][genomic_allele]['FEAT']['PL'] = pl
            # OUTPUT_LINES[genomic_pos][genomic_allele]['FEAT']['GT'] = GT_PL
            OUTPUT_LINES[genomic_pos][genomic_allele]['FEAT']['GQ'] = GQ 
            OUTPUT_LINES[genomic_pos][genomic_allele]['FEAT']['GT_str'] = GT_str

        for genomic_allele in list(OUTPUT_LINES[genomic_pos]):
            if 'GT_str' not in OUTPUT_LINES[genomic_pos][genomic_allele]['FEAT']:
                OUTPUT_LINES[genomic_pos].pop(genomic_allele)

            # print("alle", genomic_pos, genomic_allele_list, pl, Q, list(row_sums))

             



    


