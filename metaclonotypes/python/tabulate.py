"""
Script for Durable expansion of TCR-δ meta-clonotypes after BCG revaccination in humans
James, C. et al. 

Script Authored: January 13, 2021 (kmayerbl)
Reviewed: April 30, 2021 (kmayerb/)

This script searches for abundance of 2021_01_13_delta_BCG_metaclones.tsv in a 
series of bulk files.
"""

# In this file we take quantify metaclonotypes in bulk files:
import pandas as pd
import os
import pandas as pd 
from tcrdist.adpt_funcs import import_adaptive_file
from tcrdist.repertoire import TCRrep
from tcrdist.tabulate import tabulate
import math
import time
import pandas as pd
from tcrdist.repertoire import TCRrep
from collections import defaultdict 
import os 
import pandas as pd
import re

def get_safe_chunk(search_clones, bulk_clones,target = 10**7):
	"""
	This function help pick a chunk size that prevents excessive memory use,
	With two CPU, 10*7 should keep total overall memory demand below 1GB
	"""
	ideal_divisor = (search_clones * bulk_clones) / target
	if ideal_divisor < 1:
		ideal_chunk_size = search_clones
		print(ideal_chunk_size)
	else:
		ideal_chunk_size = math.ceil((search_clones)/ ideal_divisor)
		print(ideal_chunk_size)
	return ideal_chunk_size

def do_search2(file , df_search, dest, tag, path):
	
	sample_name = file.replace('.tcrdist.tsv','')
	tic = time.perf_counter()
	
	# <tr_search> tcrdist.repertoire.TCRrep object for computing distances
	tr_search = TCRrep(cell_df = df_search,
					organism = 'human',
					chains = ['delta'],
					db_file = 'alphabeta_gammadelta_db.tsv',
					compute_distances = False)
	# set cpus according to parameter above
	tr_search.cpus = 1
	df_bulk   = pd.read_csv(os.path.join(path, file), sep = '\t').rename(columns = {'cdr3_b_aa': 'cdr3_d_aa'})
	print(df_bulk)
	df_bulk   = df_bulk[['cdr3_d_aa','v_d_gene','j_d_gene','templates','productive_frequency']].rename(columns = {'templates':'count'})

	tr_bulk = TCRrep(	cell_df = df_bulk,                 
						organism = 'human', 
						chains = ['delta'], 
						db_file = 'alphabeta_gammadelta_db.tsv',
						compute_distances = False)
	
	#lines_per_file.append(tr_bulk.clone_df.shape[0]) 
	
	search_clones = tr_search.clone_df.shape[0]
	bulk_clones   = tr_bulk.clone_df.shape[0]
	# To avoid memory pressure on the system we set a target that tcrdist doesn't do more than 10M comparisons per process
	ideal_chunk_size = get_safe_chunk(tr_search.clone_df.shape[0], tr_bulk.clone_df.shape[0],target = 10**7)
	tr_search.compute_sparse_rect_distances(df = tr_search.clone_df, df2 = tr_bulk.clone_df, chunk_size = ideal_chunk_size) #(5)
	r1 = tabulate(clone_df1 = tr_search.clone_df, clone_df2 = tr_bulk.clone_df, pwmat = tr_search.rw_delta, 
		cdr3_name = 'cdr3_d_aa', v_gene_name = 'v_d_gene', j_gene_name = 'j_d_gene')

	outfile = os.path.join(dest, f"{sample_name}.{tag}.bulk_tabulation.tsv")
	print(f"WRITING: {outfile}")
	r1.to_csv(outfile, sep = '\t', index = False)
	toc = time.perf_counter()
	print(f"TABULATED IN {toc - tic:0.4f} seconds")
	del(tr_search)
	del(tr_bulk)
	#return(r1)
	return(f"{toc - tic:0.4f}s")


if __name__ == "__main__":
	# <path> where files reside
	path = os.path.join( 
		'/Volumes/T7/Seshadri', 
		'sampleExport_2020-07-22_23-33-32', 
		'tcrdist3ready')
	# <dest> where files go
	dest = os.path.join( 
	'/Volumes/T7/Seshadri', 
	'sampleExport_2020-07-22_23-33-32', 
	'tcrdist3ready', 'tabTRD')

	# <all_files> list of all files
	all_files = [f for f in os.listdir(path) if f.endswith('TRD.tcrdist.tsv')]
	search_file = '/Volumes/T7/Seshadri/2021_01_13_delta_BCG_metaclones.tsv'
	df_search = pd.read_csv(search_file, sep = '\t')
	do_search2(file = all_files[0], df_search = df_search, dest = ".", tag = "test", path = path)
	import parmap
	parmap.map(do_search2, all_files, path = path, df_search = df_search, dest = dest, tag = "tab1TRD", pm_pbar = True, pm_processes = 4)





# path = os.path.join( 
# 	'/Volumes/T7/Seshadri', 
# 	'sampleExport_2020-07-22_23-33-32', 
# 	'tcrdist3ready')
# # <all_files> list of all files
# all_files = [f for f in os.listdir(path) if f.endswith('TRD.tcrdist.tsv')]
# df_bulk = pd.read_csv(os.path.join(path, all_files[0]), sep = "\t").rename(columns = {'cdr3_b_aa':'cdr3_d_aa'})

# df_bulk = df_bulk[['cdr3_d_aa',
#                     'v_d_gene',
#                     'j_d_gene',
#                     'templates',
#                     'productive_frequency',
#                     'valid_cdr3']].\
#         rename(columns = {'templates':'count'})

# df_bulk = df_bulk[(df_bulk['v_d_gene'].notna()) & (df_bulk['j_d_gene'].notna())].reset_index()
# tr_bulk = TCRrep(cell_df = df_bulk,
#                  organism = 'human',
#                  chains = ['delta'],
#                  db_file = 'alphabeta_gammadelta_db.tsv',
#                  compute_distances = False)

# search_file = '/Volumes/T7/Seshadri/delta_BCG_metaclones.tsv'
# df_search = pd.read_csv(search_file, sep = '\t')
# df_search = df_search[['cdr3_d_aa','v_d_gene','j_d_gene','radius','regex']]
# tr_search = TCRrep(cell_df = df_search,
#                    organism = 'human',
#                    chains = ['delta'],
#                    db_file = 'alphabeta_gammadelta_db.tsv',
#                    compute_distances = False)
# tr_search.cpus = 1
# tic = time.perf_counter()
# tr_search.compute_sparse_rect_distances(df = tr_search.clone_df, df2 = tr_bulk.clone_df, chunk_size = 50, radius = 50) 
# results = tabulate(
# 	clone_df1 = tr_search.clone_df, 
# 	clone_df2 = tr_bulk.clone_df, 
# 	pwmat = tr_search.rw_delta, 
# 	cdr3_name = 'cdr3_d_aa', v_gene_name = 'v_d_gene', j_gene_name = 'j_d_gene')
# toc = time.perf_counter()
# print(f"TABULATED IN {toc - tic:0.4f} seconds")


# do_search2(file , df_search, dest, tag) 







