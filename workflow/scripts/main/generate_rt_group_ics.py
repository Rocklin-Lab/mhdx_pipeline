import os
import sys
import pandas as pd

sys.path.append(os.getcwd()+"/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx

library_info=pd.read_csv(snakemake.input[0])
name = snakemake.wildcards['name']

tp_inputs = [[tensor for tensor in snakemake.input if snakemake.config[tp][0].split('/')[-1].split('.')[0].split('_')[-1] in tensor] for tp in snakemake.config['timepoints']]

t1 = hx.TensorGenerator(name, library_info, tp_inputs, snakemake.config['timepoints'])
t1.generate_tensors()

hx.limit_write(t1.all_tp_clusters, snakemake.output[0])
