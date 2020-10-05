import os
import sys
import pandas as pd

sys.path.append(os.getcwd()+"/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx


def weak_pareto_dom_filter(po):
    out = []
    for tp in po.all_tp_clusters:
        
        tp_buffer = []
        center_dict = {}
        
        #make dict for tp int mz bins
        for i in range(len(tp[0].baseline_integrated_mz)):
            center_dict[i] = []
        
        #add ic to bin list closest to center
        for ic in tp:
            center_dict[np.round(ic.baseline_integrated_mz_com)].append(ic)
        
        #score all ics in each int_mz bin, keep only those that are best in 
        for i in range(len(tp[0].baseline_integrated_mz)): 
            int_mz_buffer = []
            score_df = pd.DataFrame().from_dict({
            "idx": [j for j in range(len(center_dict[i]))],
            "rt_ground_err": [ic.rt_ground_err for ic in center_dict[i]],
            "dt_ground_err": [ic.dt_ground_err for ic in center_dict[i]],
            "peak_err": [ic.baseline_peak_error for ic in center_dict[i]]
            })
            
            if len(score_df) > 0:
                for key in ["rt_ground_err", "dt_ground_err", "peak_err"]:
                    win_idx = int(score_df.sort_values(key).iloc[0]['idx'])
                    if win_idx not in int_mz_buffer:
                        tp_buffer.append(center_dict[i][win_idx])
                        int_mz_buffer.append(win_idx)
                
        out.append(tp_buffer)
    
    return out


library_info=pd.read_csv(snakemake.input[0])
name = snakemake.wildcards['name']
atc = hx.limit_read(snakemake.input[1])
fatc = []
for tp in atc:
    fatc.append([ic for ic in tp if ic.baseline_peak_error/ic.baseline_auc < 0.25])

p1 = hx.PathOptimizer(name, fatc, library_info, timepoints = snakemake.config['timepoints'], n_undeut_runs = len(snakemake.config[0]), old_data_dir = snakemake.config['old_data_dir'])

p1.prefiltered_ics = weak_pareto_dom_filter(p1)
p1.generate_sample_paths()

p1.optimize_paths()

p1.bokeh_plot(snakemake.output[0])

#save winner, runners, undeut_grounds, winner_scores, and rtdt com cvs
hx.limit_write(p1.winner, snakemake.output[1])
hx.limit_write(p1.runners, snakemake.output[2])
hx.limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products], snakemake.output[3])
hx.limit_write(p1.winner_scores, snakemake.output[4])
hx.limit_write([p1.rt_com_cv, p1.dt_com_cv], snakemake.output[5])



