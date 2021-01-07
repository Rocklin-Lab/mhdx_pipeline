import os
import sys
import pandas as pd

sys.path.append(os.getcwd()+"/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx

library_info = pd.read_csv(snakemake.input[0])
my_fn = snakemake.input[1]
for tp in snakemake.config['timepoints']:
	for fn in snakemake.config[tp]:
		if fn in my_fn:
			my_tp = tp

tg = hx.TensorGenerator(my_fn, my_tp, library_info)
all_ics = []
for factor in tg.DataTensor.factors:
	for ic in factor.isotope_clusters:
		all_ics.append(ic)

hx.limit_write(all_ics, snakemake.output[0])
