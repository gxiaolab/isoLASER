#!/usr/bin/python
import networkx as nx
import numpy as np
import itertools
import copy
from collections import defaultdict
from time import time

from . import gqv_software_management 

def construct_de_Bruijn_graph(PartialReads, options):
	G = nx.DiGraph()
	k = options.kmer

	for read_name, attr in PartialReads.items():
		ReadSeq = attr['seq']
		for a_j in range(len(ReadSeq) - k + 1):
			kmer_or = ReadSeq[a_j : a_j + k - 1]
			kmer_dt = ReadSeq[a_j + 1 : a_j + k]
			rname_node = (read_name, a_j)
			try:
				G[kmer_or][kmer_dt]['reads'].append(rname_node)
				G[kmer_or][kmer_dt]['r'].add(read_name)
			except:
				G.add_edge(kmer_or, kmer_dt, reads = [rname_node], r = {read_name})

	return G



def dbg_remove_lowly_covered_edges(db, REF_Index, min_Allele_Ratio):

	edges_to_remove = set()
	for node in db.nodes:
		OutNodesCountTotal = 0.0
		for node_out in db.successors(node):
			OutNodesCountTotal += len(db[node][node_out]['r'])

		for node_out in db.successors(node):
			OutNodesCount  = len(db[node][node_out]['r'])
			EdgeCountRatio = OutNodesCount/OutNodesCountTotal
			if not any(x[0] == REF_Index for x in db[node][node_out]['reads']): 
				if EdgeCountRatio <= min_Allele_Ratio: ## or OutNodesCount <= 1:
					edges_to_remove.add((node, node_out))

		InNodesCountTotal = 0.0
		for node_in in db.predecessors(node):
			InNodesCountTotal += len(db[node_in][node]['r'])

		for node_in in db.predecessors(node):
			InNodesCount   = len(db[node_in][node]['r'])
			EdgeCountRatio = InNodesCount/InNodesCountTotal
			if not any(x[0] == REF_Index for x in db[node_in][node]['reads']): 
				if EdgeCountRatio <= min_Allele_Ratio: #or InNodesCount <= 1:
					edges_to_remove.add((node_in, node))
	
	for (a, b) in edges_to_remove:
		db.remove_edge(a, b)

	del edges_to_remove


def dbg_remove_alternative_ends(all_nodes, db, node_First, node_Last):	
	for node in all_nodes:
		next_nodes = list()
		if db.out_degree(node) == 0 and db.in_degree(node) == 0:
			db.remove_node(node)
		elif db.out_degree(node) == 1 and node in db.successors(node):
			db.remove_node(node)
		elif db.in_degree(node) == 1 and node in db.predecessors(node):
			db.remove_node(node)
		elif db.in_degree(node) == 0 and node != node_First:
			if db.out_degree(node) == 0:
				next_nodes = []
			else:
				next_nodes = db.successors(node)
			db.remove_node(node)
			dbg_remove_alternative_ends(next_nodes, db, node_First, node_Last)
		elif db.out_degree(node) == 0 and node != node_Last:
			if db.in_degree(node) == 0:
				next_nodes = []
			else:
				next_nodes = db.predecessors(node)
			db.remove_node(node) 
			dbg_remove_alternative_ends(next_nodes, db, node_First, node_Last)
	


def dbg_compress_graph(db, xdb, k, prev_node, node):
	consecutive_nodes = list()
	OV = k - 1
	while db.out_degree(node) == 1 and db.in_degree(node) == 1:
		consecutive_nodes.append(node) ; a = node
		node = db.successors(node).__next__()

	if not consecutive_nodes:
		curr_node = node
	else:
		curr_node = ''.join(n[0] for n in consecutive_nodes) + consecutive_nodes[-1][1:]
	
	if (prev_node, curr_node) not in xdb.edges:
		a, b = prev_node[-OV:], curr_node[:OV]
		n_db = db[a][b]['r']
		r_db = db[a][b]['reads']

		xdb.add_edge(prev_node, curr_node, reads = r_db, r = n_db, ris = set()) 
		if node != curr_node[-OV:]:
			n_db = db[curr_node[-OV:]][node[:OV]]['r']
			r_db = db[curr_node[-OV:]][node[:OV]]['reads']
			xdb.add_edge(curr_node, node, reads = r_db, r = n_db, ris = set()) 

		for out_node in db.successors(node):
			dbg_compress_graph(db, xdb, k, node, out_node)



def dbg_get_read_paths(xdb, RG_PATHS, k):
	for (cn, nn) in xdb.edges:	
		for (rname, kmer_index) in xdb[cn][nn]['reads']:
			index_offset = len(cn) - (k - 1)
			try:
				RG_PATHS[rname][cn].append(kmer_index - index_offset)
			except:
				RG_PATHS[rname][cn] = [kmer_index - index_offset]

			# if xdb.out_degree(nn) == 0:
			try:
				RG_PATHS[rname][nn].append(kmer_index + 1)
			except:
				RG_PATHS[rname][nn] = [kmer_index + 1]



def compress_read_paths(xdb, k, RG_PATHS, PartialReads):
	OV = k - 2

	tmp1_SetOfReadIndexes = defaultdict(set)
	tmp2_SetOfReadIndexes = defaultdict(set) 

	for RG, attr in list(PartialReads.items()):

		RG_Seq = attr['seq']
		if RG in RG_PATHS:
			path_dict = RG_PATHS.pop(RG)
		else:
			PartialReads.pop(RG)
			continue

		known_nodes = set()
		for LongNode, index_list in path_dict.items():
			for r_pos in index_list:
				known_nodes.add((r_pos, LongNode))

		known_nodes = list(known_nodes)
		known_nodes.sort(key = lambda j: j[0])
		
		
		cont_read_paths = {}

		first_node = known_nodes[0][1]
		cont_read_path_length = len(first_node)
		j_0 = 0
		j_i = 0
		for j_i in range(len(known_nodes) - 1):
			(curr_i, curr_node), (next_i, next_node) = known_nodes[j_i : j_i + 2]
			if curr_node in xdb.nodes and next_node in xdb.successors(curr_node):
				cont_read_path_length += len(next_node) - OV
			else:
				cont_read_paths[cont_read_path_length] = (j_0, j_i + 1)
				j_0 = j_i + 1
				cont_read_path_length = len(next_node)

		cont_read_paths[cont_read_path_length] = (j_0, j_i + 2) 

		prev_range = [-2e9, -2e9]
		for i, (pathLen, (Max_s, Max_e)) in enumerate(cont_read_paths.items()): 
			convert_read_pos, PATH = map(list, zip(*known_nodes[Max_s : Max_e]))
			LongSeq = PATH[0] + ''.join(p[OV:] for p in PATH[1:])

			if RG == 0:
				assert LongSeq == RG_Seq, "{}\n!=\n{}".format(LongSeq, RG_Seq)
				newRG = 0
			elif len(set(PATH)) == 1 and len(PATH) > 1:
				continue
			else:
				newRG = "{}_{}".format(RG, i)

			PartialReads[newRG] = copy.deepcopy(PartialReads[RG])
			PartialReads[newRG]['Convert_Pos'] = convert_read_pos 
			PartialReads[newRG]['seq']         = LongSeq
			
			tmp1_SetOfReadIndexes[';'.join(PATH)].add(newRG)

			if Max_s < prev_range[1] and prev_range[0] < Max_e:
				raise ValueError('Overlapping ranges') 

			prev_range = [Max_s, Max_e]

		if RG != 0:
			PartialReads.pop(RG) 

	for PATH_str in sorted(tmp1_SetOfReadIndexes.keys(), key = lambda x : len(x), reverse = True):

		RGs = tmp1_SetOfReadIndexes[PATH_str]

		PATH = PATH_str.split(";")
		
		while xdb.in_degree(PATH[0]) == 1 and xdb.predecessors(PATH[0]).__next__() not in PATH \
			and not any(special_char in xdb.predecessors(PATH[0]).__next__() for special_char in ["X", "Y"]):
			
			pred = xdb.predecessors(PATH[0]).__next__()
			PATH.insert(0, pred)
			
			for RG in RGs:
				PartialReads[RG]['Convert_Pos'].insert(0, np.nan)

		while xdb.out_degree(PATH[-1]) == 1 and xdb.successors(PATH[-1]).__next__() not in PATH \
			and not any(special_char not in xdb.successors(PATH[-1]).__next__() for special_char in ["X", "Y"]):
			
			succ = xdb.successors(PATH[-1]).__next__()
			PATH.append(succ)
			
			for RG in RGs:
				PartialReads[RG]['Convert_Pos'].append(np.nan)
				

		PATH_str = ';'.join(PATH)
		tmp2_SetOfReadIndexes[PATH_str] |= RGs

	SetOfReadIndexes = []

	for Read_Index, (PATH_str, PartialReads_Set) in enumerate(tmp2_SetOfReadIndexes.items()):
		PATH = PATH_str.split(";") 
		INDEX = [0]
		icm = 0
		LongSeq = PATH[0]
		for j in range(len(PATH) - 1):
			LongSeq += PATH[j + 1][OV:]
			icm += len(PATH[j]) - OV
			INDEX.append(icm)

		### read object 
		read_obj = Read_Index_Object(INDEX, LongSeq, k, PartialReads_Set)

		for RG in PartialReads_Set:
			PartialReads[RG]['Index'] = Read_Index

			offset = np.array(INDEX) - np.array(PartialReads[RG]['Convert_Pos'])
			offset = offset[~np.isnan(offset)][0]
			PartialReads[RG]['Convert_Pos'] = list(map(int, np.array(INDEX) - offset))

		for j in range(len(PATH) - 1):
			xdb[PATH[j]][PATH[j + 1]]['ris'].add(Read_Index)

		SetOfReadIndexes.append(read_obj)

	del tmp1_SetOfReadIndexes
	del tmp2_SetOfReadIndexes

	return SetOfReadIndexes



class Read_Index_Object:

    def __init__(self, path_index, long_seq, k, PartialReads_Set):
        self.PathIndex = path_index
        self.PathIndexExt = path_index + [len(long_seq)]
        self.LongSeq = long_seq
        self.reads = PartialReads_Set
        self.k = k
        self.gStart = None
        self.gEnd = None
        self.junctions = set([])
        self.Events = defaultdict(set)
        self.EventsCount = 0
        self.mapping = []

    def PathSeq(self, return_index = False):
        array = []
        i2 = 0
        for i in range(len(self.PathIndex) - 1):
            i1, i2 = self.PathIndex[i : i + 2]
            node = self.LongSeq[i1 : i2 + self.k - 2]
            array.append(node)

        array.append(self.LongSeq[i2:])
        if return_index:
            node_index = defaultdict(list)
            for i, node in enumerate(array):
                node_index[node].append(i)
            return (array, node_index)
        else:
            return np.array(array)

    def add_mapping(self, mapping_list):
        if mapping_list:
            self.mapping = mapping_list
            self.gStart = mapping_list[0][2]
            self.gEnd = mapping_list[-1][3]

    def add_event(self, Ref_Coords, Rid_Coords, f_Seq, d_Seq):
        if len(d_Seq) - len(f_Seq) <= 20 and (d_Seq or f_Seq):
            try:
                self.Events[Rid_Coords].add(Ref_Coords)
            except:
                self.Events[Rid_Coords] = set([Ref_Coords])
            self.EventsCount += 1

    def find_position(self, POS, END):
        rid_pos, rid_end = (None, None)
        for i, (rid_s, rid_e, ref_s, ref_e, tag) in enumerate(self.mapping):
            if tag == 'REF':
                if ref_s <= POS and POS < ref_e and rid_pos is None:
                    rid_pos = (rid_s + POS - ref_s, i)
                if ref_s < END and END <= ref_e and rid_end is None:
                    rid_end = (rid_s + END - ref_s, i)
                if ref_s <= POS and END <= ref_e:
                    rid_pos, rid_end = (rid_s, i), (rid_e, i)
                    break
            elif tag == 'ALT1' and (ref_e - ref_s) > 20:
                if ref_s < POS and END < ref_e:
                    rid_pos, rid_end = (None, None)
            elif tag == 'ALT1':
                if ref_s <= POS and POS <= ref_e:
                    rid_pos = (rid_s, i)
                if ref_s <= END and END <= ref_e:
                    rid_end = (rid_e, i)

        if rid_pos is None or rid_end is None:
            if POS == END and rid_pos is not None:
                j = rid_pos[1]
                return (rid_pos[0], rid_pos[0], self.mapping[j][4])
            elif POS == END and rid_end is not None:
                j = rid_end[1]
                return (rid_end[0], rid_end[0], self.mapping[j][4])
            else:
                return (None, None, None)
        else:
            if rid_pos[1] == rid_end[1]:
                j = rid_pos[1]
                return (rid_pos[0], rid_end[0], self.mapping[j][4])
            else:
                return (rid_pos[0], rid_end[0], 'ALT1')

    def delete_events(self):
        gqv_software_management.rm_nested_dict(self.Events)



def is_contained(path_1, path_2):
    i = -1
    overlap0 = False
    while path_2[0] in path_1[i + 1: ]:
        i = path_1.index(path_2[0], i + 1)
        n = len(path_2)
        if path_2 == path_1[i : i + n]:
            overlap0 = True 
            break
    return  overlap0



def find_Source_Target_pairs3(xdb, REF_Index, SetOfReadIndexes, options):
	OV = options.kmer - 2
	REFERENCE_PATH = list(SetOfReadIndexes[REF_Index].PathSeq())

	REF_bubble_dict = defaultdict(set)  
	SRCSNK_pairs = defaultdict(set)
	
	for node in REFERENCE_PATH:
		assert node in xdb.nodes
		if xdb.out_degree(node) > 1:
			REF_bubble_dict["src"].add(node) 
		if xdb.in_degree(node) > 1:
			REF_bubble_dict["snk"].add(node) 

	for Read_Index, r in enumerate(SetOfReadIndexes):

		if Read_Index == REF_Index:
			continue		

		READ_PATH  = list(r.PathSeq())
		bubble_dict = defaultdict(set)

		for node in READ_PATH:
			if node in REF_bubble_dict["src"]:
				bubble_dict["src"].add(node)
			if node in REF_bubble_dict["snk"]:
				bubble_dict["snk"].add(node)

		for src_seq, snk_seq in itertools.product(bubble_dict["src"], bubble_dict["snk"]):
			SRCSNK_pairs[Read_Index].add((src_seq, snk_seq))

	return SRCSNK_pairs




class seq_interval_obj():
	def __init__(self, src, start, end, dummy_var = None):
		self.src = src
		if int(end) < int(start):
			raise ValueError() 
		self.start = int(start)
		self.end   = int(end)
		self.repeat = 0

	def __len__(self):
		return self.end - self.start 

	def get_seq(self, reference_sequence):
		return reference_sequence[self.start : self.end].upper()

	def add_repeat(self):
		self.repeat += 1

	def __str__(self):
		return '{0.src}|{0.start}|{0.end}|{0.repeat}'.format(self)


def cc(path_string):
	ps = path_string.split(';')
	if len(ps) > 3:
		i = 0
		new_path_string = ''
		while i < len(ps):
			src, start = ps[i].split('|')[:2]
			while i + 1 < len(ps) and ps[i + 1].startswith(src):
				i += 1
			end = ps[i].split('|')[2]
			new_path_string += '{}|{}|{};'.format(src, start, end)
			i += 1
		return new_path_string 
	else:
		return path_string



def make_read_dbg(Read_Index, Ref_Index, SetOfReadIndexes, RC):

	EVENTS = SetOfReadIndexes[Read_Index].Events
	Rid_Long_Seq = SetOfReadIndexes[Read_Index].LongSeq
	Ref_Long_Seq = SetOfReadIndexes[Ref_Index].LongSeq
	Rid_Long_Seq_ext = "ZZ" + Rid_Long_Seq + "ZZ" 

	Graph = nx.DiGraph()

	rid_nodes = set([0 , 2, len(Rid_Long_Seq_ext) - 2, len(Rid_Long_Seq_ext)])
	events_n  = 0
	events_skipped = 0

	for Rid_iv in EVENTS:
		rid_nodes.add(Rid_iv.start + 2)
		rid_nodes.add(Rid_iv.end + 2)
		events_n += len(EVENTS[Rid_iv])

	rid_nodes = list(rid_nodes)
	rid_nodes.sort()

	for i in range(len(rid_nodes) - 2):
		i1, i2, i3 = rid_nodes[i : i + 3]
		node_src = seq_interval_obj("RID", i1, i2)
		node_snk = seq_interval_obj("RID", i2, i3) 
		Graph.add_edge(str(node_src), str(node_snk), label = "READSEQ")
		
	for Rid_iv in EVENTS:
		r_i1 = rid_nodes.index(Rid_iv.start + 2)
		r_i2 = rid_nodes.index(Rid_iv.end + 2, r_i1)

		node_src = seq_interval_obj("RID", rid_nodes[r_i1 - 1], rid_nodes[r_i1])
		node_snk = seq_interval_obj("RID", rid_nodes[r_i2], rid_nodes[r_i2 + 1])
		
		for Ref_iv in EVENTS[Rid_iv]:

			s1, s2 = Rid_iv.get_seq(Rid_Long_Seq), Ref_iv.get_seq(Ref_Long_Seq)

			if s1 == s2:
				continue
			
			while str(Ref_iv) in Graph.nodes:
				Ref_iv.add_repeat() 

			if events_n >= 250 and ((s1 in s2) or (s2 in s1)):
				events_skipped += 1
				continue 

			if len(s1)*len(s2) > (1e3):
				continue 
			
			if events_n >= 250 and abs(len(Ref_iv) - (node_snk.start - node_src.end)) > 5:
				events_skipped += 1
				continue 
			
			Graph.add_edge(str(node_src), str(Ref_iv), label = "REFSEQ")
			Graph.add_edge(str(Ref_iv), str(node_snk), label = "REFSEQ")

	source = seq_interval_obj("RID", rid_nodes[0], rid_nodes[1])
	target = seq_interval_obj("RID", rid_nodes[-2], rid_nodes[-1])

	return Graph, str(source), str(target)




def Dijkstra_find_end_to_end_paths(Read_Index, Ref_Index, SetOfReadIndexes, Graph, 
		source, target, edges, RC_i):

	Ref_Long_Seq = SetOfReadIndexes[Ref_Index].LongSeq
	Rid_Long_Seq = SetOfReadIndexes[Read_Index].LongSeq

	path = dict()
	Rid_Long_Seq_ext = "ZZ" + Rid_Long_Seq + "ZZ"

	for node in Graph.nodes:
		node_obj = seq_interval_obj(*node.split('|'))
		if node_obj.src == "RID":
			node_obj_seq = node_obj.get_seq(Rid_Long_Seq_ext)
		else:
			node_obj_seq = node_obj.get_seq(Ref_Long_Seq)

		path[node] = {'SEQ' : node_obj_seq.strip('Z'), 'PATH' : {}}

	if not path:
		return {}

	path[source]['PATH'][source] = {'ref_pos' : None, 'path_seq' : path[source]['SEQ']}

	queue = set([(source, 0)])
	
	stack_counter = 0

	while queue:

		(u, up) = min(queue, key = lambda x: x[1])
		queue.remove((u, up))

		u_obj = seq_interval_obj(*u.split("|"))

		# ref_pos_set = set([vv['ref_pos'] for p, vv in path[u]['PATH'].items() ])
		
		ref_pos_set_counter = defaultdict(int)
		if len(path[u]['PATH']) * Graph.out_degree(u) > 2e5:
			for u_path, u_attr in list(path[u]['PATH'].items()):
				ref_pos = u_attr['ref_pos']
				if u_path.count("REF") == 0:
					continue
				elif ref_pos_set_counter[ref_pos] >= 15:
					path[u]['PATH'].pop(u_path)
				else:
					ref_pos_set_counter[ref_pos] += 1 


		for v in Graph.successors(u):

			v_obj = seq_interval_obj(*v.split("|"))
			v_obj_Seq = path[v]['SEQ']

			for u_path, u_attr in path[u]['PATH'].items():
				stack_counter += 1

				path_coord    = u_path + ';' + v 
				path_sequence = u_attr['path_seq'] + v_obj_Seq
				ref_pos       = u_attr['ref_pos']

				if v_obj.src == "REF" and Ref_Long_Seq[:v_obj.end].endswith(path_sequence):

					tmp_ref_pos = v_obj.end - len(path_sequence)

					if ref_pos is not None and tmp_ref_pos != ref_pos:
						continue
			
					w = Graph.successors(v).__next__()
					
					w_obj = seq_interval_obj(*w.split("|"))
					w_obj_Seq = path[w]['SEQ']
					
					assert Graph.out_degree(v) == 1 and Graph.in_degree(v) == 1
					assert w_obj.src == "RID" and u_obj.src == "RID"

					path_coord_ext    = u_path + ';' + v + ';' + w
					path_sequence_ext = u_attr['path_seq'] + v_obj_Seq + w_obj_Seq

					if Ref_Long_Seq[tmp_ref_pos:].startswith(path_sequence_ext):
						if path_coord_ext.count('REF') > 15:
							continue
						path[w]['PATH'][path_coord_ext] = {'ref_pos': tmp_ref_pos, 'path_seq': path_sequence_ext}
						edges[(v, w)].append('REF')
						queue.add((w, w_obj.start))
						
				elif v_obj.src == "RID":
					if ref_pos is None or Ref_Long_Seq[ref_pos:].startswith(path_sequence):
						path[v]['PATH'][path_coord] = {'ref_pos': ref_pos, 'path_seq': path_sequence}
						queue.add((v, v_obj.start))

		if u != target:
			path.pop(u)

		if stack_counter > int(2e6):
			gqv_software_management.print_time_stamp("WARNING: Too many recursions for this RC.\n")
			queue = set([])
	 
	for _path, _attr in list(path[target]['PATH'].items()):
		pos, seq = _attr['ref_pos'], _attr['path_seq']
		if pos is None:
			if "REF" not in _path and Rid_Long_Seq in Ref_Long_Seq:
				pos = Ref_Long_Seq.index(Rid_Long_Seq)
				path[target]['PATH'][_path]['ref_pos'] = pos
			else:
				path[target]['PATH'].pop(_path)
		else:
			assert Ref_Long_Seq[pos : pos + len(seq)] == seq

	start_to_end_paths = dict((p, path[target]['PATH'][p]['ref_pos']) for p in path[target]['PATH'])
	
	gqv_software_management.rm_nested_dict(path) 

	return start_to_end_paths


def _dbg_compress_x(nd, k):
	if len(nd) > 2*k:
		return nd[:k] + '_' + nd[-k:]
	else:
		return nd


def draw_rid_graph(xdb, graph_file_name):
	out = open("{}.txt".format(graph_file_name), "w")
	out.write('digraph G {\n')
	out.write('size="7,7"\nrankdir=TB;\n')

	for nd in xdb.nodes:
		nd = nd.replace("|","_")
		if nd.startswith("RID"):
			color = "red"
		else:
			color = "blue"
		string1 = f'\t{nd} [color={color},label={nd}];\n'
		out.write(string1)

	for nd_or, nd_dt in xdb.edges:
		nd_or = nd_or.replace("|","_")
		nd_dt = nd_dt.replace("|","_")
		label = ''
		if nd_or.startswith("RID") and nd_dt.startswith("RID"):
			color = "red"
		else:
			color = "black"

		string = f'\t{nd_or} -> {nd_dt} [color={color},label="{label}"];\n'
		out.write(string)

	out.write('}')
	out.close()


def draw_de_Bruijn_graph(xdb, k, graph_file_name, Ref_Index = 0):
	out = open("{}.txt".format(graph_file_name), "w")
	out.write('digraph G {\n')
	out.write('size="7,7"\nrankdir=TB;\n')

	for nd in xdb.nodes:
		nd_str  = _dbg_compress_x(nd, k)
		string1 = '\t{0} [label={0}_{1}];\n'.format(nd_str, len(nd) - k + 2)
		out.write(string1 )

	for nd_or, nd_dt in xdb.edges:
		N = len(xdb[nd_or][nd_dt]['r'])

		label = "{}".format(N)
		color = "blue"

		for x in xdb[nd_or][nd_dt]['reads']:
			if type(x) is tuple and x[0] == Ref_Index :
				label = "{}_REF".format(N)
				color = "red"
				break
			elif type(x) is int and x == Ref_Index:
				label = "{}_REF".format(N)
				color = "red"
				break

		# if len(xdb[nd_or][nd_dt]['reads']) == 0:
		# 	continue

		nd_or = _dbg_compress_x(nd_or, k)
		nd_dt = _dbg_compress_x(nd_dt, k)

		string = '\t{0} -> {1} [color={3},label="{2}"];\n'.format(nd_or, nd_dt, label, color)
		out.write(string)

	out.write('}')
	out.close()




