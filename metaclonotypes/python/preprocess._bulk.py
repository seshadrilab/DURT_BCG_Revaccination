#preprocess.py

# 2021-01-11 
# This is simply file conversion. We take adaptive output and grab relevant fields.
# We also subset out TRDV | TRDJ receptors from the TCRa files

import pandas as pd
import os
from progress.bar import IncrementalBar
from tcrdist.adpt_funcs import import_adaptive_file, adaptive_to_imgt

def safe_startswith(x,s):
	if isinstance(x, str):
		return x.startswith(s)
	else:
		return False
				

r = '/Volumes/T7/Seshadri/sampleExport_2020-07-22_23-33-32'
dest = os.path.join(r,"tcrdist3ready")
if not os.path.isdir(dest):
	os.mkdir(dest)

fs = [f for f in os.listdir(r) if f.endswith('.tsv')]
cache = list()
bar = IncrementalBar('Processing', max=len(fs))
for f in fs:
	f = os.path.join(r,f)
	assert os.path.isfile(f)
	if f.find("TCRa") !=-1:
		d = import_adaptive_file(f, organism = "human", chain = "alpha")
		d['subject'] == os.path.basename(f)
		d['calc_sum_temp'] =  d['templates'].sum()
		cache.append(
			{"filename":os.path.basename(f),
			'calc_sum_templates': d['templates'].sum(),
			'nrows': d.shape[0]})
		outname = f"{os.path.basename(f)}.tcrdist.tsv"
		fp_outname = os.path.join(dest, outname )
		print(f"WRITING {outname} to {fp_outname}")
		d.to_csv(fp_outname, sep = "\t", index = False)
		# subset only the deltas

		ind1 = d['v_a_gene'].apply(lambda x: safe_startswith(x = x , s = "TRDV"))
		ind2 = d['j_a_gene'].apply(lambda x: safe_startswith(x = x, s = "TRDJ"))
		d_delta = d[(ind1)|(ind2)].reset_index(drop = True)
		d_delta = d_delta.rename(columns = {'v_a_gene':'v_d_gene','j_a_gene':'j_d_gene','cdr3_a_aa':'cdr3_d_aa' })
		d_delta['calc_sum_temp'] = d_delta['templates'].sum()
		outname = f"{os.path.basename(f)}.TRD.tcrdist.tsv"
		fp_outname = os.path.join(dest, outname )
		print(f"WRITING {outname} to {fp_outname}")
		d_delta.to_csv(fp_outname, sep = "\t", index = False)

	elif f.find("TCRB") !=-1:
		d = import_adaptive_file(f, organism = "human", chain = "beta")
	else:
		raise ValueError(f"unexpected chain for {f}")
	d['subject'] == os.path.basename(f)
	d['calc_sum_temp'] =  d['templates'].sum()
	cache.append(
		{"filename":os.path.basename(f),
		'calc_sum_templates': d['templates'].sum(),
		'nrows': d.shape[0]})
	outname = f"{os.path.basename(f)}.tcrdist.tsv"
	fp_outname = os.path.join(dest, outname )
	print(f"WRITING {outname} to {fp_outname}")
	d.to_csv(fp_outname, sep = "\t", index = False)
	bar.next()
bar.finish()