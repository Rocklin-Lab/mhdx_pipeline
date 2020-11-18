import os
import sys
import copy
import math
import numpy as np
import pandas as pd

sys.path.append(os.getcwd()+"/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx

library_info=pd.read_csv(snakemake.input[0])
name = snakemake.wildcards['name']
atc = hx.limit_read(snakemake.input[1])

p1 = hx.PathOptimizer(name, atc, library_info, timepoints = snakemake.config['timepoints'], n_undeut_runs = len(snakemake.config[0]), old_data_dir = snakemake.config['old_data_dir'])
p1.optimize_paths()
p1.bokeh_plot(snakemake.output[0])

#save winner, runners, undeut_grounds, winner_scores, and rtdt com cvs
hx.limit_write(p1.winner, snakemake.output[1])
hx.limit_write(p1.runners, snakemake.output[2])
hx.limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products], snakemake.output[3])
hx.limit_write(p1.winner_scores, snakemake.output[4])
hx.limit_write([p1.rt_com_cv, p1.dt_com_cv], snakemake.output[5])



