#!/usr/bin/python
import numpy as np
import sys

FEATURE_LIST = ['INDEL_LEN', 'AB', "REP_COUNT", "ndHAP", "logAC"] 
REF_dummies  = ['REF_base_A', 'REF_base_C', 'REF_base_G', 'REF_base_T', 'REF_base_N', 'REF_base_-'] 
ALT_dummies  = ['ALT_base_A', 'ALT_base_C', 'ALT_base_G', 'ALT_base_T', 'ALT_base_N', 'ALT_base_-'] 
GT_dummies   = ['GT_1/1', 'GT_0/1', 'GT_0/0']


def modify_alleles_str(REF, ALT):
	r = [0.0 for x in REF_dummies]
	a = [0.0 for x in ALT_dummies]

	if REF and ALT and REF[0] == ALT[0]:
		REF, ALT = REF[1:], ALT[1:]
	
	if len(ALT) > 1 and len(set(ALT)) == 1:
		ALT = ALT[0]
	elif len(ALT) > 1 and len(set(ALT)) > 1:
		ALT = "N"	
	elif len(ALT) == 0:
		ALT = "-"

	if len(REF) > 1 and len(set(REF)) == 1:
		REF = REF[0]
	elif len(REF) > 1 and len(set(REF)) > 1:
		REF = "N"
	elif len(REF) == 0:
		REF = "-"

	if REF not in ['-', 'A', 'C', 'G', 'T', 'N']:
		REF = 'N'
		
	if ALT not in ['-', 'A', 'C', 'G', 'T', 'N']:
		ALT = 'N'
	
	r_i = REF_dummies.index(f'REF_base_{REF}')
	a_i = ALT_dummies.index(f'ALT_base_{ALT}')

	r[r_i] = 1.0 ; a[a_i] = 1.0 

	assert sum(a) == 1.0 and sum(r) == 1.0
	return (r , a)


def modify_gt_str(gt_string):
	a = [0.0 for x in GT_dummies]
	if np.nansum(gt_string) == 2:
		new_gt_string= '1/1'
	elif np.nansum(gt_string) == 1:
		if np.sum(np.isnan(gt_string)) == 1:
			new_gt_string = "1/1"
		else:
			new_gt_string = '0/1'
	elif np.nansum(gt_string) == 0:
		new_gt_string = '0/0'

	a_i = GT_dummies.index(f"GT_{new_gt_string}")

	a[a_i] = 1.0
	assert sum(a) == 1.0
	return a


def update_vars(CLF_dict, ALL_VARS_FEAT, options, pseudo_prob = 1e-10):
	for POS, ALT_dict in ALL_VARS_FEAT.items():
		for ALLELE, feat_dict in ALT_dict.items():
			(REF, ALT, END) = ALLELE
			ATTR = feat_dict["FEAT"]

			var_type = ATTR['variant_type']
			
			if var_type in ("INTRONIC_PART", "INTRON", "EditingSite", "REF_BLOCK"):
				ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["QUAL"] = np.nan	
			
			elif var_type == 'COMPLEX':
				if ATTR['INDEL_LEN'] > 0:
					var_type = "INSERTION"
				elif ATTR['INDEL_LEN'] < 0:
					var_type = "DELETION"
				else:
					var_type = "SNV"	
			
			else:
				FEATURE_VECTOR = [ATTR[key] for key in FEATURE_LIST] 
				
				_ref_dummies, _alt_dummies = modify_alleles_str(REF, ALT)	

				_gt_dummies = modify_gt_str(ATTR["GT"])
				
				FEATURE_VECTOR = FEATURE_VECTOR + _ref_dummies + _alt_dummies + _gt_dummies 

				
				# missing variant types ? 	
				clf = CLF_dict[var_type]["clf"]
				scl = CLF_dict[var_type]["scaler"]

				
				FEATURE_VECTOR = scl.transform([FEATURE_VECTOR])[0]
				FEATURE_VECTOR = np.array(FEATURE_VECTOR)

				idx = np.isnan(FEATURE_VECTOR) 
				FEATURE_VECTOR[idx] = scl.mean_[idx]

				try:
					p_False, p_True = clf.predict_proba([FEATURE_VECTOR])[0]
					QUAL = p_True #-10 * np.log10(p_False + pseudo_prob)
				except:
					print("problem feature vector", POS, ALLELE, FEATURE_LIST, FEATURE_VECTOR)
					QUAL = 0.0 
				
				if QUAL > 0.5 and ATTR["DP"] >= options.minCoverage and ATTR["AC"] >= options.minCountVar:
					FILTER = 'PASS'
				else:
					FILTER = 'LowQual'

				ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["QUAL"]   = QUAL
				ALL_VARS_FEAT[POS][ALLELE]["FEAT"]["FILTER"] = FILTER








