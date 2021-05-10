"""
Script for Durable expansion of TCR-δ meta-clonotypes after BCG revaccination in humans
James, C. et al. 

Script Authored: January 13, 2021 (kmayerbl)
Reviewed: April 30, 2021 (kmayerb/)

This script take a file of expanded clones and 
forms public meta-clonotypes, with fixed TCRdist radius 18. 
"""
import pandas as pd
import numpy as np
import os
from tcrdist.repertoire import TCRrep
from tcrsampler.sampler import TCRsampler
from tcrdist.adpt_funcs import import_adaptive_file, adaptive_to_imgt
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt

# Declare the project folder
project_folder = '/Volumes/T7/Seshadri/BCG'
# <fp> is full file path to data recieved Nov 11, 2020 from C.James
fp = os.path.join(project_folder, 'data-raw', '20201111_expanded_clones_postBCG.csv')
# <raw> is DataFrame of raw data 
raw = pd.read_csv(fp)
# <selections> is a dictionary that will be used to rename and select important columns
selections = {'PTID':'subject', 
'sampleA_abundance':'count_a',
'sampleB_abundance':'count_b',
'sampleA_frequency':'freq_a',
'sampleB_frequency':'freq_b',
'aminoAcid': 'cdr3_x_aa',
'vGeneName': 'v_x_gene',
'jGeneName': 'j_x_gene',
'TCRchain' : 'chain',
'TimePoint': 'timepoint'}
# <alpha_or_delta_df> temporary DataFrame
alpha_or_delta_df = raw[list(selections.keys())].rename(columns =selections)

#<alpha_df> contains only alpha seuqences
alpha_df = alpha_or_delta_df.query("chain =='alpha'").\
	copy().\
	reset_index(drop = True).\
	rename(columns = {	'cdr3_x_aa':'cdr3_a_aa',
						'v_x_gene' :'v_a_gene',
						'j_x_gene' :'j_a_gene'})
#<delta_df> contains only alpha seuqences
delta_df = alpha_or_delta_df.query("chain =='delta'").\
	copy().\
	reset_index(drop = True).\
	rename(columns = {	'cdr3_x_aa':'cdr3_d_aa',
						'v_x_gene' :'v_d_gene',
						'j_x_gene' :'j_d_gene'})

# Check fo missing v_gene_names
any_missing_delta = Counter(delta_df.v_d_gene.apply(lambda x : (x,adaptive_to_imgt['human'].get(x))))
[(k,v) for k,v in any_missing_delta.items() if k[1] is None or k[1] == "NA"]

any_missing_alpha = Counter(alpha_df.v_a_gene.apply(lambda x : (x,adaptive_to_imgt['human'].get(x))))
[(k,v) for k,v in any_missing_alpha.items() if k[1] is None or k[1] == "NA"]

# Substitute Adaptive gene names for IMGT names
delta_df.v_d_gene = delta_df.v_d_gene.apply(lambda x : adaptive_to_imgt['human'].get(x))
delta_df.j_d_gene = delta_df.j_d_gene.apply(lambda x : adaptive_to_imgt['human'].get(x))

alpha_df.v_a_gene = alpha_df.v_a_gene.apply(lambda x : adaptive_to_imgt['human'].get(x))
alpha_df.j_a_gene = alpha_df.j_a_gene.apply(lambda x : adaptive_to_imgt['human'].get(x))

alpha_df.to_csv(os.path.join(project_folder, "data"," alpha_df.tsv"), sep = "\t", index = False)
delta_df.to_csv(os.path.join(project_folder, "data"," delta_df.tsv"), sep = "\t", index = False)

# See more details here: https://tcrdist3.readthedocs.io/en/latest/#tcr-distancing
trd = TCRrep(	cell_df = delta_df, 
				organism = 'human', 
				chains = ['delta'], 
				db_file = 'alphabeta_gammadelta_db.tsv')
# Matrix of delta-chain pairwise distances
trd.pw_delta

# The tcrdist delta-chain matrix is available here and can be easily visualized:
gd = sns.clustermap(data= trd.pw_delta,
                   row_cluster=True,
                   col_cluster=True,
                   yticklabels=False,
                   xticklabels=False,
                  )

# FIND METACLONOTYPES
from tcrdist.public import _neighbors_fixed_radius
tcrsampler_delta = TCRsampler(default_background = 'ravens_human_delta_t.sampler.tsv')
trd.clone_df['radius']   = 18
trd.clone_df['neighbors'] = _neighbors_fixed_radius(pwmat = trd.pw_delta, radius = 18)
trd.clone_df['K_neighbors'] = trd.clone_df['neighbors'].apply(lambda x : len(x))
trd.clone_df['nsubject']   = trd.clone_df['neighbors'].\
    apply(lambda x: trd.clone_df['subject'].iloc[x].nunique())
trd.clone_df['qpublic']   = trd.clone_df['nsubject'].\
    apply(lambda x: x > 1)

from tcrdist.public import make_motif_logo
from tcrdist.public import _quasi_public_meta_clonotypes
qpublic_mcs = _quasi_public_meta_clonotypes(clone_df = trd.clone_df, 
											pwmat = trd.pw_delta,
											tcrsampler = tcrsampler_delta, 
											cdr3_name = 'cdr3_d_aa',
											v_gene_name = 'v_d_gene',
											nr_filter = True,
											output_html_name = "TRD_quasi_public_clones.html",
											sort_columns = ['nsubject','K_neighbors'],
											sort_ascending = False)

# Add regex for purposes of tabulation
from tcrdist.regex import _index_to_regex_str
qpublic_mcs['quasi_public_df']['regex'] = [_index_to_regex_str(ind = i, 
	clone_df = trd.clone_df.copy(), 
	pwmat = trd.pw_delta.copy(),
	col = 'cdr3_d_aa',
	max_ambiguity = 5) for i in qpublic_mcs['quasi_public_df']['neighbors']]
qpublic_mcs['quasi_public_df'][['nsubject','v_d_gene', 'j_d_gene','cdr3_d_aa','radius','regex']].\
	to_csv('2021_01_13_delta_BCG_metaclones.tsv',sep = "\t", index = True)

# For Tabulation of these feature against bulk data
# do python project/tabulation.py

