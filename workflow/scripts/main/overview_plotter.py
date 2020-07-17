import glob
import sys
import os
import zlib
import copy
import numpy as np
import pandas as pd
import seaborn as sns
import _pickle as cpickle
import matplotlib as mpl
mpl.use("AGG")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.append(os.getcwd()+"/scripts/auxiliary")
import LC_IM_MS_TensorAnalysis as ta

from bokeh.models import HoverTool, ColorBar, Text, Div
from bokeh.io import push_notebook, show, output_notebook, output_file, save
from bokeh.models.callbacks import CustomJS
from bokeh.models.sources import ColumnDataSource, CDSView
from bokeh.models.filters import Filter, GroupFilter, IndexFilter
from bokeh.plotting import figure
from bokeh.layouts import gridplot, column
from bokeh.transform import linear_cmap
from bokeh.palettes import Spectral6

import ipdb

#Make and display mass-added plot
def overview_mass_added_plotter(source, tooltips):
    p = figure(title='Mean Added Mass over HDX Timeseries, Colored by RTxDT Error in ms ', plot_height=350, plot_width=1300, background_fill_color='whitesmoke', tools='pan,wheel_zoom,box_zoom,hover,save,reset,help', tooltips=tooltips)
    err_low, err_high = min(source.data['rtxdt_rmse']), max(source.data['rtxdt_rmse'])
    mapper = linear_cmap(field_name='rtxdt_rmse', palette=Spectral6, low=err_low, high=2)
    color_bar = ColorBar(color_mapper=mapper['transform'], width=10,  location=(0,0))
    p.multi_line(xs='timepoints', ys='mz_coms', source=source, color=mapper, hover_color='red', line_width=1.5)
    p.xaxis.axis_label = "Timepoint Index"
    p.yaxis.axis_label = "Added-Mass Units"
    p.add_layout(color_bar, 'right')
    return p

def overview_error_plotter(source, tooltips):
    p = figure(title='RMSEs of RT-Group RT and DT Center-of-Mass Errors', plot_height=350, plot_width=1300, background_fill_color='whitesmoke', tools='pan,wheel_zoom,box_zoom,hover,save,reset,help', tooltips=tooltips)
    p.circle(x="rt_rmse", y="dt_rmse", size = 5, alpha = 0.4, hover_color="red", hover_alpha=1, source=source)
    p.xaxis.axis_label = "RMSE of RT-Error (ms)"
    p.yaxis.axis_label = "RMSE of DT-Error (ms)"
    return p



#Input handling depends on naming conventions
winner_paths = sorted([path for path in snakemake.input if "winner.cpickle.zlib" in path])
undeut_ground_paths = sorted([path for path in snakemake.input if "undeut_grounds" in path])
winner_score_paths = sorted([path for path in snakemake.input if "winner_scores" in path])
rtdt_cv_paths = sorted([path for path in snakemake.input if "rtdt_com_cvs" in path])

pnames = sorted(set(["_".join(fn.split("/")[-1].split('_')[:3]) for fn in winner_paths])) 
snames = sorted(set(["_".join(fn.split("/")[-1].split('_')[:4]) for fn in winner_paths]))

#Making dataset for overview plots

out = {}
for winner_path in winner_paths:
    #Create values needed for bokeh portion of overview plot:
        #Array of timepoint values in order (range(len(winner)))
        #Array of baseline_int_mz_com values in order
        #Mean of (rt bin err * (0.07) x dt abs err)
        #RT-cluster mean RT value in Sname
    
    pname = "_".join(winner_path.split("/")[-1].split("_")[:3])
    sname = "_".join(winner_path.split("/")[-1].split("_")[:4])
    rt_cluster = winner_path.split("/")[-1].split("_")[3]
    winners = ta.limit_read(winner_path)
    timepoints = list(range(len(winners)))

    
    mz_coms = []
    rt_errs = []
    dt_errs = []
    rtxdt_errs = []
    for ic in winners:
        mz_coms.append(ic.baseline_integrated_mz_com)
        rt_errs.append((ic.rt_ground_err*0.07)**2)
        dt_errs.append(ic.dt_ground_err**2) #dt error is in ms
        rtxdt_errs.append(((ic.rt_ground_err*0.07)*ic.dt_ground_err)**2)#0.07 approximately scales rt-bins to ms

    rt_rmse = np.sqrt(np.mean(rt_errs))
    dt_rmse = np.sqrt(np.mean(dt_errs))
    rtxdt_rmse = np.sqrt(np.mean(rtxdt_errs))
    
    
    if pname in out.keys():
        out[pname].append((pname, sname, rt_cluster, timepoints, mz_coms, rtxdt_rmse, rt_rmse, dt_rmse))
   
    else:
        out[pname] = []
        out[pname].append((pname, sname, rt_cluster, timepoints, mz_coms, rtxdt_rmse, rt_rmse, dt_rmse))

data = []
[data.append(species_tup) for key in out.keys() for species_tup in out[key]]

source_frame = pd.DataFrame(data, columns= ["pname", "sname", "mean_rt", "timepoints", "mz_coms", "rtxdt_rmse", "rt_rmse", "dt_rmse"])

source = ColumnDataSource(source_frame)

tooltips = [
    ("Protein-Name", "@pname"),
    ("Mean-RT of RT-Group", "@mean_rt"),
    ("RT COM RMSE to Undeuterated", "@rt_rmse"),
    ("DT COM RMSE to Undeuterated", "@dt_rmse")
]

### Collect data for distplots that can't be done in bokeh ###
score_keys = [
    "int_mz_std_rmse", 
    "baseline_peak_error", 
    "delta_mz_rate", 
    "dt_ground_rmse", 
    "dt_ground_fit", 
    "rt_ground_fit", 
    "rt_ground_rmse", 
    "auc_ground_rmse"
    ]

new_scores = dict.fromkeys(score_keys)
score_weights = dict.fromkeys(score_keys)

for key in new_scores.keys():
    new_scores[key] = []

new_scores_switch = None
for path in winner_score_paths:
    score_dict = ta.limit_read(path)
    if new_scores_switch is None:
        for key in score_weights.keys():
            score_weights[key] = score_dict[key][0]
        new_scores_switch = True
    for key in new_scores.keys():
        new_scores[key].append(score_dict[key][1])

#keep fit closest to 1 for all charges in rt-group
idotps = [] #idotp = isotopic dot product: the dot product of the observed and theoretical added mass distributions of an isotopic cluster
for path in undeut_ground_paths:
    fit_dict =  ta.limit_read(path)[1]
    best_fit = 0
    for key in fit_dict.keys():
        if fit_dict[key] > best_fit:
            best_fit = fit_dict[key]
    idotps.append(best_fit)

#Add to score dict
new_scores['idotp'] = idotps

#Add rtdt com cvs to new_scores
new_scores['rt_com_cvs'], new_scores['dt_com_cvs']  = [], []
for path in rtdt_cv_paths:
    tup = ta.limit_read(path)
    new_scores['rt_com_cvs'].append(tup[0])
    new_scores['dt_com_cvs'].append(tup[1])

new_score_frame = pd.DataFrame.from_dict(new_scores)
new_score_frame['type'] = ['New' for i in range(len(new_score_frame))]

score_frame = copy.copy(new_score_frame)

#Make old-data dist-plot datasources if called for in config
if snakemake.config['old_data_dir'] is not None:
    
    #output individual stat values for displotting as 1D arrays
    old_scores = dict.fromkeys(['idotp', 'delta_mz_rate'])
    for key in old_scores.keys():
        old_scores[key] = [] 

    #get old files and populate dict of lists
    for file in glob.iglob("data/old_data*.pickle"): #TODO: This needs to look for files in data/old_data now
        with pickle.loads(open(file, 'rb').read()) as ts:
            old_scores['idotp'].append(ts['fit_to_theo_dist'])
            old_scores['delta_mz_rate'].append(ts['delta_mz_rate'])

    old_scores['type'] = ['Old' for i in range(len(old_scores))]

    #make mixed df for new and old data
    mixed_scores = []
    for key in old_scores.keys():
        mixed_scores.append(pd.Series(old_scores[key], name=key))
    for key in new_scores.keys():
        mixed_scores.append(pd.Series(new_scores[key], name=key))
    score_frame = pd.concat(mixed_scores, axis=1)
    
    g = sns.PairGrid(new_score_frame, hue="type")
    g.map_diag(sns.distplot)
    g.map_upper(plt.scatter)
    g.map_lower(sns.kdeplot)
    g.add_legend();
    plt.savefig("plots/"+snakemake.config['run_name']+"_overview_pairplot.png")
else:
    g = sns.PairGrid(score_frame)
    g.map_diag(sns.distplot)
    g.map_upper(plt.scatter)
    g.map_lower(sns.kdeplot)
    plt.savefig("plots/"+snakemake.config['run_name']+"_overview_pairplot.png")

#Pairplot can now be accessed as a source in a bokeh Div(<img>) obj

final = column(Div(text="<h1 style='margin-left: 250px;'>Overview of "+snakemake.config['run_name']+"</h1>"), overview_mass_added_plotter(source, tooltips), overview_error_plotter(source, tooltips), Div(text="<img source='plots/"+snakemake.config['run_name']+"_overview_pairplot.png' width='1300' height='1300>'"))

output_file(snakemake.output[0], mode='inline')

save(final)




