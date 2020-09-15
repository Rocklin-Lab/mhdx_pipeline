import os
import sys
import glob
import zlib
import math
import copy
import pickle
import pymzml
import importlib.util
import numpy as np
import pandas as pd
import seaborn as sns
import _pickle as cpickle
import matplotlib.pyplot as plt
from importlib import reload
from collections import defaultdict as ddict
from scipy.ndimage.filters import gaussian_filter
from tensorly.decomposition import non_negative_parafac as nnp

sys.path.append(os.getcwd()+"/workflow/scripts/auxiliary/")
import hxtools
import LC_IM_MS_TensorAnalysis as hx

ins = copy.copy(snakemake.input)

library_info = pd.read_csv(ins.pop(0))
i = int(ins[0].split('/')[-1].split('_')[0])

hd_mass_diff = 1.006277
c13_mass_diff = 1.00335
low_mass_margin = 10
high_mass_margin = 17
ppm_radius = 30

#make hxprot representation per sequence
name = library_info.iloc[i]['name']
charge = library_info.iloc[i]['charge']
seq = library_info.loc[library_info['name'] == name]['sequence'].values[0]
hx_prot = hxtools.hxprot(seq=seq)
dot_products = []

#Compute DataTensor init inputs
max_peak_center = len(seq)
total_isotopes = max_peak_center+high_mass_margin 
total_mass_window = low_mass_margin+total_isotopes
est_peak_gaps = [0] + list(np.linspace(c13_mass_diff, hd_mass_diff, 7)) + [hd_mass_diff for x in range(total_isotopes - 8)]
cum_peak_gaps = np.cumsum(est_peak_gaps)

mz_lows = library_info['obs_mz'].values[i] - (low_mass_margin/library_info['charge'].values[i])
mz_highs = library_info['obs_mz'].values[i] - (total_isotopes/library_info['charge'].values[i])
mz_centers = library_info['obs_mz'].values[i] + (cum_peak_gaps/library_info['charge'].values[i])
low_lims = mz_centers * ((1000000.0 - ppm_radius)/1000000.0)
high_lims = mz_centers * ((1000000.0 + ppm_radius)/1000000.0)

#find rt_group undeut tensors
undeut_fns = [fn for fn in glob.glob("resources/tensors/"+str(i)+"_*") if "UN" in fn]

print("Debug: Start undeut factorizations")
count = 1
undeut_ics = []
for fn in undeut_fns:
    #read mzml 
    print("Loop: "+str(count)+", Debug 1")
    output = hx.limit_read(fn)
    print("Loop: "+str(count)+", Debug 2")
    #make DataTensor
    newDataTensor = hx.DataTensor(
                            source_file = fn,
                            tensor_idx = i,
                            timepoint_idx = 0, 
                            name = name, 
                            total_mass_window = total_mass_window,
                            n_concatenated = 1,
                            charge_states = [library_info["charge"].values[i]], 
                            rts = output[0], 
                            dts = output[1], 
                            seq_out = output[2], 
                            int_seq_out = None
                            )
    print("Loop: "+str(count)+", Debug 3")
    print("Tensor Dims: "+str(np.shape(newDataTensor.full_grid_out)))
    
    newDataTensor.lows = np.searchsorted(newDataTensor.mz_labels, low_lims)
    newDataTensor.highs = np.searchsorted(newDataTensor.mz_labels, high_lims)
    newDataTensor.factorize(gauss_params=(3,1))
    print("Loop: "+str(count)+", Debug 4")
    #Hacky, clean
    [undeut_ics.append(ic) for i in range(len(newDataTensor.factors)) for ic in newDataTensor.factors[i].isotope_clusters]

    print("Loop: "+str(count)+", Debug 5")
    count += 1

for ic in undeut_ics:
    df = pd.DataFrame(ic.baseline_integrated_mz, columns = ['major_species_integrated_intensities'])
    fit = hx_prot.calculate_isotope_dist(df)
    ic.undeut_ground_dot_product = fit
    dot_products.append((fit, ic.charge_states))

if len(dot_products) > 0:
	best_idotp = max(np.asarray(dot_products)[:,0])
else:
	best_idotp = 0

print("Debug: Loop Complete")

hx.limit_write(pd.DataFrame({'index': [i], 'name': [name], 'charge': [charge], 'idotp': [best_idotp]}), snakemake.output[0])