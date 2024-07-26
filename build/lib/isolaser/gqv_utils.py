#!/usr/bin/python

import sys
import numpy as np
import scipy.sparse
import copy
import itertools
import datetime
import random
from collections import defaultdict, Counter
from time import time
import pyfaidx 
import HTSeq
import networkx as nx
import warnings

from . import gqv_software_management as SOFTMAN
from . import gqv_dbg_utils
from . import gqv_glm
from . import gqv_vartool


warnings.filterwarnings('ignore')


def Process_de_Bruijn_Graph(PartialReads, options, db):
    k = options.kmer

    node_First  = PartialReads[0]['seq'][:k - 1]
    node_Last   = PartialReads[0]['seq'][-(k - 1):]
    node_Second = PartialReads[0]['seq'][1 : k]
    
    all_nodes = list(db.nodes)[:]
    
    ## essential given the high number of errors
    gqv_dbg_utils.dbg_remove_lowly_covered_edges(db, 0, options.maxSeqErrorRate)

    gqv_dbg_utils.dbg_remove_alternative_ends(all_nodes, db, node_First, node_Last) 

    xdb = nx.DiGraph()

    gqv_dbg_utils.dbg_compress_graph(db, xdb, k, node_First, node_Second) 

    xdb.name = db.name 

    RG_PATHS = defaultdict(dict)

    # gqv_dbg_utils.draw_de_Bruijn_graph(xdb, k, xdb.name, Ref_Index = 0) ; print('drawn')

    gqv_dbg_utils.dbg_get_read_paths(xdb, RG_PATHS, k)
    
    SetOfReadIndexes = gqv_dbg_utils.compress_read_paths(xdb, k, RG_PATHS,
                                                            PartialReads)

    SOFTMAN.rm_nested_dict(RG_PATHS)

    return (xdb, SetOfReadIndexes)




def assign_vars_to_reads(src_snk_pairs, Ref_Index, SetOfReadIndexes, PartialReads, 
        k, REFPOS):
    
    memoization_vars = defaultdict(set)
    memoization_trims = dict()
    SetOfVariants = [[]]

    map_rate = [[0, 0, 0], [0, 0, 0]]
    times = []

    for ronda in (1, 2):
        for Read_Index in range(len(SetOfReadIndexes)):
            
            if SetOfReadIndexes[Read_Index].mapping:
                continue            

            n = len(SetOfReadIndexes[Read_Index].reads)

            t00 = time() 
            mapping2 = [] 
            for tmp_var_list in SetOfVariants:
                mapping2 = _align_read_with_vars_2(Read_Index, SetOfReadIndexes, 
                                REFPOS, tmp_var_list)
                if mapping2: 
                    map_rate[ronda - 1][0] += 1
                    SetOfReadIndexes[Read_Index].add_mapping(mapping2)
                    break

            t0 = time() 
            
            if SetOfReadIndexes[Read_Index].mapping:
                continue
            
            try:
                bubble_list_tmp = src_snk_pairs.pop(Read_Index)
            except KeyError:
                bubble_list_tmp = {}

            # SOFTMAN.print_time_stamp(f"\t\t1. debug {Read_Index}  bubbles={len(bubble_list_tmp)}")
            
            _match_src_snk_to_reads(Read_Index, Ref_Index, bubble_list_tmp, SetOfReadIndexes, 
                                    memoization_trims, memoization_vars, ronda = ronda, REFPOS = REFPOS)
        
            t1 = time()
            # SOFTMAN.print_time_stamp(f"\t\t2. debug {Read_Index} evnets={len(SetOfReadIndexes[Read_Index].Events)}")

            Graph, source, target = gqv_dbg_utils.make_read_dbg(Read_Index, Ref_Index, 
                                                                SetOfReadIndexes, 1)
            # SOFTMAN.print_time_stamp(f"\t\t3. debug {Read_Index} {len(SetOfReadIndexes)}")
            SetOfReadIndexes[Read_Index].delete_events()
            edges = defaultdict(list)
            t2 = time()
            
            all_paths = gqv_dbg_utils.Dijkstra_find_end_to_end_paths(Read_Index, Ref_Index, 
                                                                    SetOfReadIndexes, Graph, 
                                                                    source, target, edges, 1)
            # SOFTMAN.print_time_stamp(f"\t\t4. debug {Read_Index} paths={len(all_paths)}")                                              
            t3 = time()
            RI_alignment, var_list = _find_best_read_path(Read_Index, Ref_Index, all_paths, 
                                                            REFPOS, k, SetOfReadIndexes, 
                                                            PartialReads)

            if not RI_alignment and Read_Index != Ref_Index:
                map_rate[ronda - 1][2] += 1
            else:
                map_rate[ronda - 1][1] += 1

            if var_list:
                SetOfVariants.append(var_list)

            SetOfReadIndexes[Read_Index].add_mapping(RI_alignment)
            t4 = time()

            times.append([t4-t3, t3-t2, t2-t1, t1-t0, t0-t00]) 

                    

    SOFTMAN.rm_nested_dict(memoization_vars)
    SOFTMAN.rm_nested_dict(memoization_trims)
    SetOfVariants = list(itertools.chain.from_iterable(SetOfVariants))
    return SetOfVariants



def _src_snk_comb(src_i_list, snk_i_list, src_seq, snk_seq, k, common_nodes):

    effective_src_i_list = src_i_list[:]
    effective_snk_i_list = snk_i_list[:]

    if len(src_i_list) > 5:
        for i in range(len(src_i_list) - 1):
            if src_i_list[i] - src_i_list[i - 1] == 1 and src_i_list[i + 1] - src_i_list[i] == 1:
                effective_src_i_list[i] = 'NA'
        effective_src_i_list = [x for x in effective_src_i_list if x != 'NA']

    if len(snk_i_list) > 5:
        for i in range(1, len(snk_i_list) - 1):
            if snk_i_list[i] - snk_i_list[i - 1] == 1 and snk_i_list[i + 1] - snk_i_list[i] == 1:
                effective_snk_i_list[i] = 'NA'
        effective_snk_i_list = [x for x in effective_snk_i_list if x != 'NA']


    if 'N' in src_seq or 'N' in snk_seq:
        combs = []

    elif len(effective_snk_i_list) > 5 and len(effective_src_i_list) > 5:
        space1 = [effective_snk_i_list[i + 1] - effective_snk_i_list[i] for i in range(len(effective_snk_i_list) - 1)]
        space2 = [effective_src_i_list[i + 1] - effective_src_i_list[i] for i in range(len(effective_src_i_list) - 1)]

        if abs(np.mean(space1) - np.mean(space2)) < 1:
            combs = []
            for i in range(len(effective_snk_i_list[:-1])):
                snk_i = effective_snk_i_list[i + 1]
                if snk_i >= effective_src_i_list[0]:
                    combs.append((effective_src_i_list[0], snk_i))
        else:
            combs = []
            for s_i, s_j in itertools.product(effective_src_i_list, effective_snk_i_list):
                delta = s_j - s_i
                if delta >= 0 and delta < 5*k and not any(s_j >= s and s >= s_i for s in common_nodes):
                    combs.append((s_i, s_j))
    else:
        combs = []
        for s_i, s_j in itertools.product(effective_src_i_list, effective_snk_i_list):
            if s_j >= s_i and not any(s_j >= s and s >= s_i for s in common_nodes):
                combs.append((s_i, s_j))

    return combs



def _match_src_snk_to_reads(Read_Index, REF_Index, bubble_list_tmp, SetOfReadIndexes, 
        memoization_trims, memoization_vars, ronda, REFPOS = None):

    k = SetOfReadIndexes[REF_Index].k
    OV = k - 2
    Ref_Seq_path, Ref_node_index = SetOfReadIndexes[REF_Index].PathSeq(return_index=True)
    Ref_Seq_index = SetOfReadIndexes[REF_Index].PathIndex
    Ref_Long_Seq  = SetOfReadIndexes[REF_Index].LongSeq
    
    Rid_Seq_path, Rid_node_index = SetOfReadIndexes[Read_Index].PathSeq(return_index=True)
    Rid_Seq_index = SetOfReadIndexes[Read_Index].PathIndex
    Rid_Long_Seq  = SetOfReadIndexes[Read_Index].LongSeq
    
    memoization_read_haplotypes = defaultdict(list)

    if Rid_Long_Seq in Ref_Long_Seq:
        bubble_list_tmp = {}
        return None
    

    common_nodes = set(Ref_Seq_path) & set(Rid_Seq_path)
    common_nodes = [x for x in common_nodes if len(x) > (k - 1)]

    common_nodes_ref = []
    common_nodes_rid = [] 

    for (src_seq, snk_seq) in bubble_list_tmp:
        assert len(src_seq) == k - 1 and len(snk_seq) == k - 1
        ref_src_i = Ref_node_index[src_seq]
        ref_snk_i = Ref_node_index[snk_seq]
        rid_src_i = Rid_node_index[src_seq]
        rid_snk_i = Rid_node_index[snk_seq]
        ref_combs = _src_snk_comb(ref_src_i, ref_snk_i, src_seq, snk_seq, k, common_nodes_ref)
        rid_combs = _src_snk_comb(rid_src_i, rid_snk_i, src_seq, snk_seq, k, common_nodes_rid)

        # SOFTMAN.print_time_stamp(f"\t\t\t\t\t src,snk={(src_seq, snk_seq)} n= {ref_combs} {rid_combs}")

        for ref_i, ref_j in ref_combs:
            f1 = Ref_Seq_index[ref_i] + OV
            f2 = Ref_Seq_index[ref_j] + OV
            SrcSnk_Ref_Seq = Ref_Long_Seq[f1 : f2]

            for rid_i, rid_j in rid_combs:

                r1 = Rid_Seq_index[rid_i] + OV
                r2 = Rid_Seq_index[rid_j] + OV
                SrcSnk_Rid_Seq = Rid_Long_Seq[r1 : r2]              

                hap_1 = Ref_Long_Seq[:f1] + SrcSnk_Rid_Seq + Ref_Long_Seq[f2:]
                hap_2 = Rid_Long_Seq[:r1] + SrcSnk_Ref_Seq + Rid_Long_Seq[r2:]

                if hap_2 in memoization_read_haplotypes[hap_1]:
                    continue
                else:
                    memoization_read_haplotypes[hap_1].append(hap_2)

                try:
                    V_t = memoization_trims[SrcSnk_Rid_Seq, SrcSnk_Ref_Seq]
                except:
                    V_t = gqv_vartool.trim_sequences(SrcSnk_Rid_Seq, SrcSnk_Ref_Seq, first = 'sufix')
                    memoization_trims[SrcSnk_Rid_Seq, SrcSnk_Ref_Seq] = V_t

                Pfx, Sfx, d_Seq, f_Seq = V_t
 
                if abs(len(d_Seq) - len(f_Seq)) > 5*k:
                    continue

                memoization_vars[src_seq, snk_seq, SrcSnk_Ref_Seq, f1, f2].add(V_t)

                RdP = gqv_dbg_utils.seq_interval_obj('RID', r1 + Pfx, r2 - Sfx)
                RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)
                
                SetOfReadIndexes[Read_Index].add_event(RfP, RdP, f_Seq, d_Seq)
                
                memoization_vars[src_seq, snk_seq, SrcSnk_Ref_Seq, f1, f2].add(V_t)
                
    for memo_key in memoization_vars:

        src_seq, snk_seq, SrcSnk_Ref_Seq, f1, f2 = memo_key
        ref_src_i = Ref_node_index[src_seq]
        ref_snk_i = Ref_node_index[snk_seq]
        rid_src_i = Rid_node_index[src_seq]
        rid_snk_i = Rid_node_index[snk_seq]

        vars_from_comb = memoization_vars[memo_key] 
        for rid_i in rid_src_i:
            if rid_i >= len(Rid_Seq_index) - 1:
                continue
            r1_tmp = Rid_Seq_index[rid_i] + OV
            for Pfx, Sfx, d_Seq, f_Seq in vars_from_comb:
                r1 = r1_tmp + Pfx
                r2 = r1 + len(d_Seq)
                if r2 <= len(Rid_Long_Seq) and d_Seq == Rid_Long_Seq[r1 : r2] and r2 + Sfx >= len(Rid_Long_Seq) - 1:
                    
                    hap_1 = Ref_Long_Seq[:f1 + Pfx] + d_Seq + Ref_Long_Seq[f2 - Sfx:]
                    hap_2 = Rid_Long_Seq[:r1] + f_Seq + Rid_Long_Seq[r2:]
                    if hap_2 in memoization_read_haplotypes[hap_1]:
                        continue
                    else:
                        memoization_read_haplotypes[hap_1].append(hap_2)
                    RdP = gqv_dbg_utils.seq_interval_obj('RID', r1, r2)
                    RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)
                    

        for rid_j in rid_snk_i:
            if rid_j <= 0:
                continue
            r2_tmp = Rid_Seq_index[rid_j] + OV
            for Pfx, Sfx, d_Seq, f_Seq in vars_from_comb:
                r2 = r2_tmp - Sfx
                r1 = r2 - len(d_Seq)
                if r1 >= 0 and d_Seq == Rid_Long_Seq[r1:r2] and r1 - Pfx - OV <= 1:
                    hap_1 = Ref_Long_Seq[:f1 + Pfx] + d_Seq + Ref_Long_Seq[f2 - Sfx:]
                    hap_2 = Rid_Long_Seq[:r1] + f_Seq + Rid_Long_Seq[r2:]
                    if hap_2 in memoization_read_haplotypes[hap_1]:
                        continue
                    else:
                        memoization_read_haplotypes[hap_1].append(hap_2)
                    RdP = gqv_dbg_utils.seq_interval_obj('RID', r1, r2)
                    RfP = gqv_dbg_utils.seq_interval_obj('REF', f1 + Pfx, f2 - Sfx)




def _find_best_read_path(Read_Index, Ref_Index, paths, REFPOS, k, SetOfReadIndexes, PartialReads):  
    min_editd_1 = sys.maxsize
    min_editd_2 = 200.0 
    final_mapping  = list()
    final_var_list = list()

    memoization_alignments = dict()
 
    # Junctions    = SetOfReadIndexes[Read_Index].junctions
    Rid_Long_Seq = SetOfReadIndexes[Read_Index].LongSeq
    Ref_Long_Seq = SetOfReadIndexes[Ref_Index].LongSeq
    Rid_Long_Seq_ext = 'ZZ' + Rid_Long_Seq + 'ZZ'
     
    path_keys = list(paths.keys())
    path_keys.sort(key = lambda p: p.count('REF'))

    dijkstra_mem = defaultdict(lambda : sys.maxsize)

    for i, path_coord_tmp in enumerate(path_keys):        
        path_editd_1 = path_coord_tmp.count('REF') * 0.1
        path_editd_2 = path_editd_1
        path_vars    = []
        path_n_vars  = 0

        bad_path = False
        path_coord = path_coord_tmp.split(';') 

        for j, tmp_node in enumerate(path_coord):
            if j == 0 or j == len(path_coord) - 1:
                continue

            prev_node = path_coord[j - 1]
            next_node = path_coord[j + 1]
            
            if tmp_node.startswith('REF'):
                Rid_Src = gqv_dbg_utils.seq_interval_obj(*prev_node.split('|'))
                Rid_Snk = gqv_dbg_utils.seq_interval_obj(*next_node.split('|'))
                Ref_Var = gqv_dbg_utils.seq_interval_obj(*tmp_node.split('|'))
                Ref_Seq = Ref_Var.get_seq(Ref_Long_Seq)
                Rid_Seq = Rid_Long_Seq_ext[Rid_Src.end : Rid_Snk.start]

                # print("decomposing..", Ref_Seq, Rid_Seq, Ref_Var.start, Ref_Var.end)
                RPos, REnd = Rid_Src.end - 2, Rid_Snk.start - 2
                
                VAR = gqv_vartool.convert_to_var(Ref_Seq, Rid_Seq, Ref_Var.start, Ref_Var.end, 
                                                 0, REnd - RPos, Read_Index, REFPOS)
                
                if str(VAR) in memoization_alignments:
                    VAROBJ_ED = memoization_alignments[str(VAR)]
                else:
                    VAROBJ_ED = gqv_vartool.calculate_edit_distance(VAR, REFPOS)
                    memoization_alignments[str(VAR)] = VAROBJ_ED  

                VAROBJ_ED_copy = copy.deepcopy(VAROBJ_ED)

                for _var in VAROBJ_ED_copy[2]:
                    _var.shift_Rpos(RPos)
         
                path_editd_1 += VAROBJ_ED_copy[0]
                path_editd_2 += VAROBJ_ED_copy[1]
                path_vars    += VAROBJ_ED_copy[2]
                path_n_vars  += VAROBJ_ED_copy[3]

                if path_editd_2 < dijkstra_mem[(REnd, VAR.g_end)]:
                    dijkstra_mem[(REnd, VAR.g_end)] = path_editd_2 
                else:
                    bad_path = True
                    break 

                if path_editd_2 > min_editd_2:
                    bad_path = True
                    break 
               
        if bad_path:
            continue

        elif path_editd_2 <= min_editd_2: 
            reference_position = paths[path_coord_tmp]

            path_vars = [v for v in path_vars if v.get_variant_type() != "NO_VARIANT"]

            mapping = _align_read_with_vars(Rid_Long_Seq, REFPOS, reference_position, path_vars)

            # better_overlap, min_map_overlap = _mapping_overlap(SetOfReadIndexes, Read_Index, Ref_Index,
            #                                                     mapping, PartialReads, min_map_overlap)
            # if better_overlap:
            min_editd_1 = path_editd_1
            min_editd_2 = path_editd_2
            final_var_list = path_vars
            final_mapping = mapping

    return (final_mapping, final_var_list)



def _align_read_with_vars_2(Read_Index, SetOfReadIndexes, REFPOS, final_var_list):

    mapping = [] 
    rid_seq = SetOfReadIndexes[Read_Index].LongSeq
    ref_seq = REFPOS.Seq_ext 

    rpos_list = []
    for var in final_var_list:
        rpos = var.g_start - REFPOS.gi.start 
        assert ref_seq[rpos : rpos + var.L_ref()] == var.REF, "{} != {}".format(ref_seq[rpos : rpos + var.L_ref()], var.REF)
        rpos_list.append(rpos)
    
    rpos_list.sort()
    final_var_list.sort(key = lambda v : v.g_start)
    for j, rpos in enumerate(rpos_list[::-1]):
        var = final_var_list[-1 * (j + 1)]
        ref_seq = ref_seq[:rpos] + var.ALT + ref_seq[rpos + var.L_ref():] 
    
    if rid_seq in ref_seq:
        ref_start = ref_seq.index(rid_seq)
        rid_start = 0
        for var, ref_end in zip(final_var_list, rpos_list):
            match_len = ref_end - ref_start

            if ref_end < ref_start:
                continue 
            if rid_start + match_len > len(rid_seq):
                break 

            _m_1 = [rid_start, 
                    rid_start + match_len, 
                    ref_start, 
                    ref_end, "REF"]
            _m_2 = [rid_start + match_len, 
                    rid_start + match_len + var.L_alt(), 
                    ref_end, 
                    ref_end + var.L_ref(), "ALT1"]

            mapping.extend([_m_1, _m_2])
            var.ReadPos[Read_Index] = tuple(_m_2[: 2])

            rid_start = rid_start + match_len + var.L_alt()
            ref_start = ref_end + var.L_ref()
              
        mapping.append([rid_start, len(rid_seq), ref_start, ref_start + len(rid_seq) - rid_start, "REF"])
        
        for i in range(len(mapping)):
            mapping[i][2] = REFPOS.genome_pos(mapping[i][2])
            mapping[i][3] = REFPOS.genome_pos(mapping[i][3])

    return mapping 



def _align_read_with_vars(Rid_Long_Seq, REFPOS, reference_position, final_indel_list):
    mapping = []
    mapped_read_pos = REFPOS.genome_pos(reference_position) 
       
    prev_REF_end = mapped_read_pos
    prev_RID_end = 0

    for variant in sorted(final_indel_list, key = lambda x : x.g_start):
        RID_pos, RID_end = list(variant.ReadPos.values())[0] 
        RID_block = Rid_Long_Seq[prev_RID_end : RID_pos].strip("Z")
        REF_block = REFPOS.get_sequence(prev_REF_end, variant.g_start) 
        assert RID_block == REF_block, "{} != {} {}".format(RID_block, REF_block, str(variant))
        mapping.append((prev_RID_end, RID_pos, prev_REF_end, variant.g_start, 'REF'))
        mapping.append((RID_pos, RID_end, variant.g_start, variant.g_end, 'ALT1'))
        prev_RID_end = RID_end
        prev_REF_end = variant.g_end

    RID_block = Rid_Long_Seq[prev_RID_end:].strip("Z")
    REF_block = REFPOS.get_sequence(prev_REF_end, prev_REF_end + len(RID_block))
    assert RID_block == REF_block, '{} != {}'.format(len(RID_block), len(REF_block))
    mapping.append((prev_RID_end, len(Rid_Long_Seq), 
                    prev_REF_end, prev_REF_end + len(RID_block), 'REF'))

    return mapping

    

def overlap_vars_and_reads(SetOfVariants, SetOfReadIndexes, PartialReads, REFPOS):

    erc_blocks = HTSeq.GenomicArrayOfSets([REFPOS.gi.chrom], stranded = False)
    erc_blocks[REFPOS.gi] += "RC"
    
    SetOfVariants.sort(key = lambda x : x.g_start)

    for variant in SetOfVariants:

        erc_blocks[variant.make_gi()] += "VAR" #str(variant)

        for Read_Index in range(len(SetOfReadIndexes)):

            if Read_Index in variant.ReadPos:
                Allele = 'ALT'
                RI_pos, RI_end = variant.ReadPos[Read_Index]
            else:
                RI_pos, RI_end, Allele = SetOfReadIndexes[Read_Index].find_position(variant.rep_region.start, 
                                                                                    variant.rep_region.end)
            if RI_pos is None or RI_end is None:
                continue

            _a0, _b0 = (None, None)
            for i, ri_i in enumerate(SetOfReadIndexes[Read_Index].PathIndex):
                if RI_pos >= ri_i:
                    _a0, _a1 = i, RI_pos - ri_i
                if RI_end >= ri_i:
                    _b0, _b1 = i, RI_end - ri_i

            if _a0 is None or _b0 is None:
                continue

            for RG in SetOfReadIndexes[Read_Index].reads:

                if RG == 0:
                    continue
            
                variant.RGPos[RG] = Allele

    ## Ref blocks as variants
    REF_BLOCKS = defaultdict(dict)  
    
    for PartialRead, PR_attr in PartialReads.items():
        PR_s, PR_e = PR_attr['geno_blocks']
        if PR_s == PR_e:
            continue
        PR_gi = HTSeq.GenomicInterval(REFPOS.gi.chrom, PR_s, PR_e)
        erc_blocks[PR_gi] += PartialRead

    for iv, var_set in erc_blocks.steps():
        if "RC" not in var_set: ## Outside of Read Cluster
            continue

        if "VAR" in var_set: ## Interval contains variant
            continue

        for RG in var_set:
            if RG not in (0, "RC"):
                REF_BLOCKS[iv][RG] = "REF" 

 
        # for Read_Index in range(len(SetOfReadIndexes)):

        #     RI_pos, RI_end, Allele = SetOfReadIndexes[Read_Index].find_position(iv.start, iv.end)
            
        #     if RI_pos is None or RI_end is None:
        #         continue

        #     _a0, _b0 = (None, None)
        #     for i, ri_i in enumerate(SetOfReadIndexes[Read_Index].PathIndex):
        #         if RI_pos >= ri_i:
        #             _a0, _a1 = i, RI_pos - ri_i
        #         if RI_end >= ri_i:
        #             _b0, _b1 = i, RI_end - ri_i

        #     if _a0 is None or _b0 is None:
        #         continue

        #     s, e = SetOfReadIndexes[Read_Index].gStart , SetOfReadIndexes[Read_Index].gEnd
        #     print("cov?", iv, Read_Index, len(SetOfReadIndexes[Read_Index].reads),(s, e), (RI_pos, RI_end, Allele), (_a0, _b0))

        #     for RG in SetOfReadIndexes[Read_Index].reads:
        #         if RG == 0:
        #             continue
        #         REF_BLOCKS[iv][RG] = Allele 
    
    del erc_blocks 
    return SetOfVariants, REF_BLOCKS 


def get_motif(ref, alt):
    motif = ref[1:] + alt[1:]
    set_motif = set(motif)

    bases = set('ACGT')

    if len(set_motif) == 1:
        return motif[0]
    elif set_motif & bases == set(['A', 'T']):
        return "AT"
    elif set_motif & bases == set(['C', 'G']):
        return "CG"  
    else:
        return 'NA'


def get_variant_attr(SetOfVariants, REFPOS, Gene_Id):

    all2int = {'REF': 0, 'ALT': 1, 'ALT1': 0}
    OUTPUT_LINES  = defaultdict(dict)
    var_skip = [0, 0]
    for variant in SetOfVariants:

        VAR_READ_AND_ALLELES = dict()
        
        while variant.RGPos:
            Read_Name, Allele  = variant.RGPos.popitem() 
        
            ALLELEINT = all2int[Allele]
     
            ReadPairName = ':'.join(Read_Name.split(':')[:-1]) 
            
            if ReadPairName in VAR_READ_AND_ALLELES:
                if VAR_READ_AND_ALLELES[ReadPairName] == ALLELEINT:
                    pass
                else:
                    VAR_READ_AND_ALLELES.pop(ReadPairName)
            else:
                VAR_READ_AND_ALLELES[ReadPairName] = ALLELEINT

        if not VAR_READ_AND_ALLELES:
            continue

        ### Allele counts
        allele_counts = Counter(VAR_READ_AND_ALLELES.values())

        AC = allele_counts[1]
        RC = allele_counts[0]
        DP = AC + RC
        AB = round(AC/DP , 3)

        if AB < 0.05:
            var_skip[1] += 1 
            continue
        else:
            var_skip[0] += 1
   
        motif = get_motif(variant.REF, variant.ALT)

        if variant.length() == 0 :
            var_len = 1
        else:
            var_len = 1/variant.length()

        Feature_Vector = {'FILTER': 'PASS',
                         'REP_COUNT': min(variant.rep_count, 20),
                         'REP_COORDS': variant.rep_region,
                         'REP_MOTIF': motif, 
            	         'variant_type' : variant.get_variant_type(),
                         'is_editing': variant.is_editing_site(REFPOS.strand), 
                         'variant_id': str(variant),
            	         'ndHAP': np.log(1e-5), 
                         'logAC': np.log(min(AC + 1, 100)),
            	         'AC': AC,
                         'RC': RC,
            	         'DP': DP,
            	         'AB': AB,
                         'GT': [np.nan, np.nan],
                         'GQ': 0,
                         'PL': ['.', '.', '.'], 
                         'QUAL': 0.5,
            	         'INDEL_LEN': var_len,
                         'n_PLOIDY' : 2,
                         'PG': 'NA', 
                         'OVER_PLOIDY': False}

        anchor_pos    = (variant.chrom, Gene_Id, variant.g_start + 1, variant.get_variant_group())
        anchor_allele = (variant.REF, variant.ALT, variant.g_end)
        OUTPUT_LINES[anchor_pos][anchor_allele] = {"READS" : VAR_READ_AND_ALLELES, "FEAT" : Feature_Vector}

    return OUTPUT_LINES


def get_reference_attr(REF_blocks, REFPOS, Gene_Id):

    OUTPUT_LINES  = defaultdict(dict)

    for iv, read_dict in REF_blocks.items():

        VAR_READ_AND_ALLELES = dict()
        
        for Read_Name, Allele in read_dict.items(): 
        
            ALLELEINT = 0
     
            ReadPairName = ':'.join(Read_Name.split(':')[:-1]) 
            
            if ReadPairName in VAR_READ_AND_ALLELES:
                if VAR_READ_AND_ALLELES[ReadPairName] == ALLELEINT:
                    pass
                else:
                    VAR_READ_AND_ALLELES.pop(ReadPairName)
            else:
                VAR_READ_AND_ALLELES[ReadPairName] = ALLELEINT

        if not VAR_READ_AND_ALLELES:
            continue

        ### Allele counts
        allele_counts = Counter(VAR_READ_AND_ALLELES.values())

        AC = allele_counts[1]
        RC = allele_counts[0]
        DP = AC + RC
        AB = round(AC/DP , 3)
        REF_base = REFPOS.get_sequence(iv.start, iv.start + 1, adapter = False)

        Feature_Vector = {'FILTER': 'REF_BLOCK',
                         'REP_COUNT': 0,
                         'REP_COORDS': iv,
                         'REP_MOTIF': 'NA', 
            	         'variant_type' : 'REF_BLOCK',
                         'is_editing': False, 
                         'variant_id': str(iv),
            	         'ndHAP': np.log(1e-5), 
                         'logAC': np.log(min(AC + 1, 100)),
                         'QUAL': 0, 
            	         'AC': AC,
                         'RC': RC,
            	         'DP': DP,
            	         'AB': AB,
                         'GT': [0, 0],
                         'GQ': 0,
                         'PL': ['.', '.', '.'], 
                         'QUAL': np.nan,
            	         'INDEL_LEN': iv.length,
                         'n_PLOIDY' : 2,
                         'PG': 'NA', 
                         'OVER_PLOIDY': False}

        anchor_pos    = (iv.chrom, Gene_Id, iv.start + 1, 'REF_BLOCK')
        anchor_allele = (REF_base, '', iv.end)
        OUTPUT_LINES[anchor_pos][anchor_allele] = {"READS" : VAR_READ_AND_ALLELES, "FEAT" : Feature_Vector}

    return OUTPUT_LINES


def read_genome_fasta(genome_fasta_file):
    fasta = pyfaidx.Faidx(genome_fasta_file)
    return fasta


def write_feat_file(var_list, FeatFile, write_header):
    if write_header:
        Feat_Handle = open(FeatFile, 'w')
    else:
        Feat_Handle = open(FeatFile, 'a')

    Title = False

    for genomic_pos in sorted(var_list):
        for genomic_allele, line in var_list[genomic_pos].items():
            (REF, ALT, END) = genomic_allele

            if not Title:
                feat_list = list(line['FEAT'].keys())
                feat_list.sort()
                if write_header:
                    Feat_Handle.write('\t'.join(feat_list + gqv_glm.REF_dummies + gqv_glm.ALT_dummies + gqv_glm.GT_dummies) + '\n')
                Title = True

            GT_string = line['FEAT']['GT']

            if abs(line['FEAT']['INDEL_LEN']) > 20:
                continue
           
            r_dummies, a_dummies = gqv_glm.modify_alleles_str(REF, ALT)
            gt_dummies = gqv_glm.modify_gt_str(GT_string) 

            outline = [ line['FEAT'][feat] for feat in feat_list ] + r_dummies + a_dummies + gt_dummies
            Feat_Handle.write('\t'.join(map(str, outline)) + '\n')

    Feat_Handle.close()




def write_readcluster_file(ALL_READ_CLUSTERS, outfile, write_header):
    if write_header:
        f = open(outfile, 'w')
        f.write("Index\tChrom\tStart\tEnd\tStrand\tMax_Coverage\n")
    else:
        f = open(outfile, 'a')

    #Read_Clusters_Genes[(i, CHROM, Gene_start, coord)] = Read_Clusters_Blocks
    for RCG, RCB_list in ALL_READ_CLUSTERS:
        RCG_index, CHROM, RCG_start, RCG_end = RCG
        for RCB in RCB_list:
            RCB_start, RCB_end = RCB
            outline = [RCG_index, CHROM, RCG_start, RCG_end, RCB_start, RCB_end]
            f.write("\t".join(map(str, outline)) + "\n")
    f.close()


def write_mutinfo_file(mi_outlines, outfile, write_header):
    # H1e;H1i;H2e;H2i
    if write_header:
        f = open(outfile, 'w')
        f.write("#chrom\tstart\tend\tstrand\tGene_name\tLinkage_Type\tPG\tAMI\tchi_p\tCov\tPSI\tAR\tallelic_counts\tlabel\tDelta_PSI\n")
    else:
        f = open(outfile, 'a')

    for mi_outline in mi_outlines:
        f.write(mi_outline + "\n")
    f.close()



def write_vcf_file(var_list, search_regions, genome_fasta, VCF_file, sample_name, options, write_header = False):
    
    if write_header:
        VCF_handle = open(VCF_file, 'w')

        genome_faidx = read_genome_fasta(genome_fasta)

        VCF_handle.write('##fileformat=VCFv4.2\n')
        VCF_handle.write('##fileDate={}\n'.format(datetime.date.today()))
        VCF_handle.write('##CL=python {}">\n'.format(' '.join(sys.argv)))
        
        for search_region in sorted(search_regions):
            chrom = search_region[0]
            chrom_len = genome_faidx.index[chrom]['rlen']
            VCF_handle.write(f'##contig=<ID={chrom},length={chrom_len}>\n')

        VCF_handle.write('##FILTER=<ID=PASS,Description="Variant passes all filters">\n')
        VCF_handle.write('##FILTER=<ID=REF_BLOCK,Description="Non variant region">\n')
        VCF_handle.write('##FILTER=<ID=LowQual,Description="Variant does not pass one or more filtering criteria">\n')
        VCF_handle.write('##INFO=<ID=varType,Number=1,Type=String,Description="Variant type">\n')
        VCF_handle.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        VCF_handle.write('##INFO=<ID=TandemRep,Number=1,Type=Integer,Description="Count of repeat sequence adjacent to the variant">\n')
        VCF_handle.write('##INFO=<ID=ndHAP,Number=1,Type=Float,Description="Normalized dHap difference">\n')
        VCF_handle.write('##INFO=<ID=GeneId,Number=1,Type=String,Description="Gene id">\n')
        VCF_handle.write('##GVCFBlock=minGQ=0(inclusive),maxGQ=5(exclusive)\n')
        VCF_handle.write('##GVCFBlock=minGQ=20(inclusive),maxGQ=60(exclusive)\n')
        VCF_handle.write('##GVCFBlock=minGQ=5(inclusive),maxGQ=20(exclusive)\n')
        VCF_handle.write('##GVCFBlock=minGQ=60(inclusive),maxGQ=2147483647(exclusive)\n')
        VCF_handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of reads overlapping the variant position">\n')
        VCF_handle.write('##FORMAT=<ID=AC,Number=R,Type=Integer,Description="Number of reads containing the variant allele">\n')
        VCF_handle.write('##FORMAT=<ID=AB,Number=R,Type=Float,Description="Allelic balance of the variant">\n')
        VCF_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus genotype">\n')
        VCF_handle.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        VCF_handle.write('##FORMAT=<ID=PG,Number=1,Type=String,Description="Phasing Group">\n')
        VCF_handle.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype Likelihoods">\n')
        VCF_handle.write('##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">\n')
        VCF_handle.write('##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">\n')
        VCF_handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(sample_name))
    
    else:
        VCF_handle = open(VCF_file, 'a')

    
    for genomic_pos in sorted(var_list, key = lambda v : v[2]):
        
        (CHROM, Gene_Id, POS, var_group) = genomic_pos

        if not var_list[genomic_pos]:
            continue 

        ### Sort alleles by descending qual

        alleles = sorted(var_list[genomic_pos], 
                         key = lambda x : var_list[genomic_pos][x]['FEAT']['QUAL'], 
                         reverse = True)

        LONG_REF = ''     
        for genomic_allele in alleles:
            (REF, ALT, END) = genomic_allele
            if len(REF) > len(LONG_REF):
                    LONG_REF = REF


        QUAL = [line_attr['FEAT']['QUAL'] for line_attr in var_list[genomic_pos].values()]
        QUAL = np.mean(QUAL) * 100

        if QUAL > 50:
            FILTER = "PASS"
        else:
            FILTER = '.'
                
        ## ALT alleles

        ALTs = []
        for (REF, ALT, END) in var_list[genomic_pos]:
            if ALT == "" and LONG_REF[len(REF): ] == "": 
                ALT = '.'
            else:
                ALT = ALT + LONG_REF[len(REF): ]
            ALTs.append(ALT)

        ALTs = ','.join(ALTs + ['<NON_REF>']) 

     

        ## ndHAP
        ndHAP = [line_attr['FEAT']['ndHAP'] for line_attr in var_list[genomic_pos].values()]
        ndHAP = sum(ndHAP)

        ## Counts
        AC = [line_attr['FEAT']['AC'] for line_attr in var_list[genomic_pos].values()]
        AC = ','.join([f'{ac:.0f}' for ac in AC])

        DP = [line_attr['FEAT']['DP'] for line_attr in var_list[genomic_pos].values()]
        DP = max(DP)

        PG = [line_attr['FEAT']['PG'] for line_attr in var_list[genomic_pos].values()]
        PG = PG[0]

        AB = [line_attr['FEAT']['AB'] for line_attr in var_list[genomic_pos].values()]
        
        PL = [line_attr['FEAT']['PL'] for line_attr in var_list[genomic_pos].values()][0]

        GT = [line_attr['FEAT']['GT_str'] for line_attr in var_list[genomic_pos].values()][0]
        GQ = [line_attr['FEAT']['GQ'] for line_attr in var_list[genomic_pos].values()][0]

        if DP < options.minCoverage: 
            continue 
  
        if len(ALTs.split(',')) == 1:
            PL = PL 
        else: 
            PL = list(PL) + [PL[0], PL[1], PL[0]]
            

        var_is_mnp = False
        for alt in ALTs.split(','):
            if len(LONG_REF) == len(alt) and len(LONG_REF) > 1:
                var_is_mnp = True

        multi_allelic = len(ALTs.split(',')) > 2

        if var_is_mnp or multi_allelic:
            print('skipping', POS, LONG_REF, ALTs)
            continue

        PL = ','.join(map(str, PL))
        AB = ','.join([f'{ab:.2f}' for ab in AB])
            
        ## Rep count
        REP_COUNT = [line_attr['FEAT']['REP_COUNT'] for line_attr in var_list[genomic_pos].values()]
        REP_COUNT = max(REP_COUNT)

        ## Writing 
        outline_tmp_basic  = f'{CHROM}\t{POS}\t.\t{LONG_REF}\t{ALTs}\t{QUAL:.2f}\t{FILTER}\t' 
        outline_tmp_info   = f'END={END};varType={var_group};TandemRep={REP_COUNT:.0f}\t'
        outline_tmp_format = f'GT:GQ:PGT:DP:MIN_DP:AC:AB:PL:PG\t{GT}:{GQ}:{GT}:{DP}:{DP}:{AC}:{AB}:{PL}:{Gene_Id}_{PG}\n'

        VCF_handle.write(outline_tmp_basic + outline_tmp_info + outline_tmp_format)

    VCF_handle.close()



