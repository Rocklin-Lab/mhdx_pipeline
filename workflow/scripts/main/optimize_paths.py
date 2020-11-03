import os
import sys
import copy
import math
import numpy as np
import pandas as pd

sys.path.append(os.getcwd()+"/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx

def atc_rel_height_rebound(atc, bound=15):
    atc_buf = []
    for tp in atc:
        tp_buf = []
        for ic in tp: 
            int_mz = ic.baseline_integrated_mz
            baseline = max(int_mz)*0.15 #TODO HARDCODE
            center = int(ic.baseline_integrated_mz_com)
            i, j = center, center
            cutoff = int_mz[center]*(bound/100)
            
            """
            plt.plot(int_mz)
            plt.axvline(center)
            plt.show()
            plt.close()
            """
            
            while (center-i <= 10 and i-1 != -1):
                i -= 1
                if int_mz[i] < cutoff:
                    break
            while (j-center <= 10 and j+1 != len(int_mz)): 
                j += 1
                if int_mz[j] < cutoff:
                    break 
                    
            if ic.concat_dt_idxs is not None:
                tp_buf.append(
                    hx.IsotopeCluster(
                        charge_states = ic.charge_states, 
                        factor_mz_data = ic.factor_mz_data, 
                        source_file = ic.source_file,
                        tensor_idx = ic.tensor_idx, 
                        timepoint_idx = ic.timepoint_idx, 
                        n_factors = ic.n_factors, 
                        factor_idx = ic.factor_idx, 
                        cluster_idx = ic.cluster_idx, 
                        low_idx = ic.lows[i], 
                        high_idx = ic.highs[j], 
                        lows = ic.lows, 
                        highs = ic.highs, 
                        grate = ic.grate, 
                        rts = ic.rts, 
                        dts = ic.dts, 
                        max_rtdt = ic.max_rtdt, 
                        outer_rtdt = ic.outer_rtdt, 
                        box_dist_avg = ic.box_dist_avg, 
                        abs_mz_low = ic.abs_mz_low, 
                        n_concatenated = ic.n_concatenated, 
                        concat_dt_idxs = ic.concat_dt_idxs,
                        total_mass_window = ic.total_mass_window
                    )
                )
                
            else:
                tp_buf.append(
                    hx.IsotopeCluster(
                        charge_states = ic.charge_states, 
                        factor_mz_data = copy.deepcopy(ic.factor_mz_data), 
                        source_file = ic.source_file,
                        tensor_idx = ic.tensor_idx, 
                        timepoint_idx = ic.timepoint_idx, 
                        n_factors = ic.n_factors, 
                        factor_idx = ic.factor_idx, 
                        cluster_idx = ic.cluster_idx, 
                        low_idx = ic.lows[i]-math.ceil(ic.box_dist_avg/2), 
                        high_idx = ic.lows[j]-math.ceil(ic.box_dist_avg/2), 
                        lows = ic.lows, 
                        highs = ic.highs, 
                        grate = ic.grate, 
                        rts = ic.rts, 
                        dts = ic.dts, 
                        max_rtdt = ic.max_rtdt, 
                        outer_rtdt = ic.outer_rtdt, 
                        box_dist_avg = ic.box_dist_avg, 
                        abs_mz_low = ic.abs_mz_low, 
                        n_concatenated = ic.n_concatenated,
                        concat_dt_idxs = None,
                        total_mass_window = ic.total_mass_window
                    )
                )
        atc_buf.append(tp_buf)
    return atc_buf

library_info=pd.read_csv(snakemake.input[0])
name = snakemake.wildcards['name']
atc = hx.limit_read(snakemake.input[1])

ratc = atc_rel_height_rebound(atc)

p1 = hx.PathOptimizer(name, ratc, library_info, timepoints = snakemake.config['timepoints'], n_undeut_runs = len(snakemake.config[0]), old_data_dir = snakemake.config['old_data_dir'])
p1.optimize_paths()
p1.bokeh_plot(snakemake.output[0])

#save winner, runners, undeut_grounds, winner_scores, and rtdt com cvs
hx.limit_write(p1.winner, snakemake.output[1])
hx.limit_write(p1.runners, snakemake.output[2])
hx.limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products], snakemake.output[3])
hx.limit_write(p1.winner_scores, snakemake.output[4])
hx.limit_write([p1.rt_com_cv, p1.dt_com_cv], snakemake.output[5])



