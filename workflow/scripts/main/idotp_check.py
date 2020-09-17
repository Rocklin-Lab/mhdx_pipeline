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

#Tensor label-generation constants
hd_mass_diff = 1.006277
c13_mass_diff = 1.00335

#config params
low_mass_margin = 10
high_mass_margin = 17
ppm_radius = 30

ins = copy.copy(snakemake.input)
library_info = pd.read_csv(ins.pop(0))
library_info['Drift Time MS1'] = library_info['im_mono'] / 200.0 * 13.781163434903
i = int(ins[0].split('/')[-1].split('_')[0])

#make hxprot representation per sequence
name = library_info.iloc[i]['name']
charge = library_info.iloc[i]['charge']
seq = library_info.loc[library_info['name'] == name]['sequence'].values[0]
hx_prot = hxtools.hxprot(seq = seq)
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
gauss_undeut_ics = []
DTs = []
charge_ic_int_mzs = []
for fn in ins:
    #read mzml 
    output = hx.limit_read(fn)

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

    newDataTensor.lows = np.searchsorted(newDataTensor.mz_labels, low_lims)
    newDataTensor.highs = np.searchsorted(newDataTensor.mz_labels, high_lims)
    #newDataTensor.factorize()

    #Hacky, clean
    #[undeut_ics.append(ic) for i in range(len(newDataTensor.factors)) for ic in newDataTensor.factors[i].isotope_clusters]

    newDataTensor.factorize(gauss_params=(3,1))
    [gauss_undeut_ics.append(ic) for i in range(len(newDataTensor.factors)) for ic in newDataTensor.factors[i].isotope_clusters]
    DTs.append(newDataTensor)

    count += 1

for ic in gauss_undeut_ics:
    df = pd.DataFrame(ic.baseline_integrated_mz, columns = ['major_species_integrated_intensities'])
    fit = hx_prot.calculate_isotope_dist(df)
    ic.undeut_ground_dot_product = fit
    dot_products.append((fit, ic.charge_states, ic.baseline_integrated_mz))
    charge_ic_int_mzs.append(ic.baseline_integrated_mz)

if len(dot_products) > 0:
    best_idotp_idx = np.argmax(np.asarray(dot_products)[:,0])
    best_idotp = max(np.asarray(dot_products)[:,0])
else:
    best_idotp = 0

factors = [factor.integrated_mz_data for tensor in DTs for factor in tensor.factors] 
factor_frames = [[factor[:i]for i in range(8,12)] for factor in factors]
best_factor_idotp = [max([hx_prot.calculate_isotope_dist(pd.DataFrame(frame, columns = ['major_species_integrated_intensities']))for factor in factor_frames for frame in factor ]), np.argmax([hx_prot.calculate_isotope_dist(pd.DataFrame(frame, columns = ['major_species_integrated_intensities']))for factor in factor_frames for frame in factor ])]
#print("Debug: Loop Complete")


#best undeut peak w/ theor undeut
plt.figure()
plt.plot(hx_prot.calculate_theo_isotope_dist())
plt.plot(dot_products[best_idotp_idx][2])
plt.show()
plt.close()

#all ics found w/ theor undeut
plt.figure()
plt.plot(hx_prot.isotope_dist)
for ic in gauss_undeut_ics:
    plt.plot(ic.baseline_integrated_mz)
plt.show()
plt.close()

#all factor int mzs w/ theor undeut
for dt in DTs:
    plt.figure()
    plt.plot(hx_prot.isotope_dist)
    for factor in dt.factors:
        plt.plot(factor.integrated_mz_data)
    plt.show()
    plt.close()

n_factors = [dt.n_factors for dt in DTs]

print("Index: "+str(i)+", IdotP: "+str(best_idotp)+", Best Factor IdotP: "+str(best_factor_idotp))    
hx.limit_write(pd.DataFrame({'index': [i], 'name': [name], 'charge': [charge], 'idotp': [best_idotp],"factor_idotp": [best_factor_idotp], "best": [dot_products[best_idotp_idx][2]], "ic_mzs": [charge_ic_int_mzs], "factor_mzs": [factors], "n_factors": [n_factors]}), snakemake.output[0])