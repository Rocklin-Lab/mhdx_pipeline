"""
Usage:
snakemake.input: <undeuterated_timepoint.isotopes> from imtbx

Reads IMTBX formatted scans and outputs a .csv of M/z, retention times, and drift times associated with known proteins in sample mixture. 
These values are used by the stg1.py script to identify relevent scans to include in a protein DataTensor.
"""
import importlib.util
hxtools_spec = importlib.util.spec_from_file_location("hxtools.py", "workflow/scripts/auxiliary/hxtools.py")
molmass_spec = importlib.util.spec_from_file_location("molmass.py", "workflow/scripts/auxiliary/molmass.py")
hxtools = importlib.util.module_from_spec(hxtools_spec)
molmass = importlib.util.module_from_spec(molmass_spec)
hxtools_spec.loader.exec_module(hxtools)
molmass_spec.loader.exec_module(molmass)

import sys
import time
import copy
import math
import Bio.PDB
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import *
from Bio.SeqUtils import ProtParam
from scipy import signal
from sklearn.cluster import DBSCAN
from sklearn.metrics import pairwise_distances
from scipy.stats import gaussian_kde

matplotlib.use('Agg')

import ipdb


### DEFINITIONS ###
    
def plotcluster(i=0):
    plt.figure(figsize=(16,3))
    plt.subplot(141)
    plt.plot(testq[clusters == i]['RT'], testq[clusters == i]['mz_mono'])
    plt.title('mz range = %.4f' % (max(testq[clusters == i]['mz_mono']) - min(testq[clusters == i]['mz_mono'])))
    
    plt.subplot(142)
    plt.plot(testq[clusters == i]['RT'], testq[clusters == i]['cluster_im'])
    plt.title('im range = %.1f' % (max(testq[clusters == i]['cluster_im']) - min(testq[clusters == i]['cluster_im'])))
   
    plt.subplot(143)
    plt.plot(testq[clusters == i]['RT'], testq[clusters == i]['ab_cluster_total'])

    plt.subplot(144)
    plt.plot(testq[clusters == i]['RT'], testq[clusters == i]['cluster_corr'])
    plt.savefig("../../plots/"+"_RT_cluster_plots.png")
    
def kde_plot(sum_df, outpath):
    mykde=gaussian_kde(sum_df['ppm'])
    xs=np.linspace(-50,50,10000)
    xs[np.argmax(mykde.evaluate(xs))]
    sns.distplot(sum_df['ppm'])
    plt.xlim([-50,50])
    plt.grid()
    plt.plot(xs, mykde(xs))
    plt.savefig(outpath)

#mix doesn't serve any purpose in the pipeline TODO
def getnear(x, charge = None, mix = None, ppm = 50):
    subdf = allseq
    if mix != None:
        subdf = allseq[allseq['mix'] == mix]
    if charge != None:
        low, high = ((x*charge)-(1.00727 * charge)) * ((1000000 - ppm)/1000000), ((x*charge)-(1.00727 * charge)) * ((1000000 + ppm)/1000000)
        mlow, mhigh = allseq['MW'] > low, allseq['MW'] < high
        tempdf = allseq[mlow & mhigh].sort_values('MW')[['MW','mix','name','len','sequence']]
        tempdf['plus%s' % int(charge)] = [(q + (1.00727 * charge) )/charge for q in tempdf['MW']]
        tempdf['ppm'] = ['%.1f' % ((1.0 - (q / x)) * 1000000) for q in tempdf['plus%s' % int(charge)] ]
        tempdf['abs_ppm'] = [np.abs(((1.0 - (q / x)) * 1000000)) for q in tempdf['plus%s' % int(charge)] ]
        return tempdf[['plus%s' % int(charge), 'ppm', 'abs_ppm', 'MW', 'mix', 'name','len','sequence']]
    else:
        low, high = x-window, x+window
        mlow, mhigh = allseq['MW'] > low, allseq['MW'] < high
        tempdf = subdf[mlow & mhigh].sort_values('MW')[['MW','mix','name','len','sequence']]
        return tempdf

def cluster_df(testq, ppm = 50, adjusted = False):
    sum_data=[]
    for c in range(0, max(testq['cluster']) + 1):
        cluster_df= testq[testq['cluster'] == c]
        charge=np.median(cluster_df['charge'])
        if adjusted:
            mz=np.average(cluster_df['mz_mono_fix_round'], weights = cluster_df['ab_cluster_total'])
        else:
            mz=np.average(cluster_df['mz_mono'], weights = cluster_df['ab_cluster_total'])
        RT=np.average(cluster_df['RT'], weights = cluster_df['ab_cluster_total'])
        im=np.average(cluster_df['im_mono'], weights = cluster_df['ab_cluster_total'])

        near = getnear(mz, charge = charge, mix = 2, ppm = ppm)

        if len(near) > 0:
            sum_data.append([near['name'].values[0],
                             RT,
                             im,
                             sum(cluster_df['ab_cluster_total']),
                             near['MW'].values[0],
                             charge,
                             near['plus%s' % int(charge)].values[0],
                             mz,
                             near['ppm'].values[0],
                             near['abs_ppm'].values[0],
                             c])
        if len(near) > 1:
            display(near)
    sum_df=pd.DataFrame(sum_data)
    sum_df.columns=['name','RT','im_mono','ab_cluster_total','MW','charge','expect_mz','obs_mz','ppm','abs_ppm','cluster']
    sum_df['ppm'] = [float(x) for x in sum_df['ppm']]
    return sum_df

def find_offset(sum_df):
    """
    Returns suspected systemic ppm error and width of poi peak of run data from sum_df. 
    Assumes protein of interest within +/- 50ppm of 0ppm. 
    Selects closest peak to 0 if sufficiently prominent.
    """
    #Maybe make this save the distplot too.
    ppm_dist = sns.distplot(sum_df['ppm'].values).get_lines()[0].get_data()
    peaks = sp.signal.find_peaks(ppm_dist[1])[0]
    xs, ys = ppm_dist[0][peaks], ppm_dist[1][peaks]
    #If lowest ppm peak is also highest frequency within window, we hope this will be the common case in our 50 ppm initial window
    try:
        xs[np.argmin(abs(xs))] == ys[np.argmax(ys)]
    except:
        ipdb.set_trace()
    if xs[np.argmin(abs(xs))] == ys[np.argmax(ys)]:
        return xs[np.argmin(abs(xs))]
        #Lowest ppm peak is not most prominent, determine relative height of lowest ppm peak 
    else: 
        #If peak closer to zero is less than half the height of the more prominent peak, check larger peak's ppm
        if ys[np.argmin(abs(xs))] < ys[np.argmax(ys)]/2:
            #If most prominent peak is heuristically close to 0, or lowest ppm peak is relatively very small (10% of major peak): use big peak
            if xs[np.argmax(ys)] < 25 or ys[np.argmin(abs(xs))] < ys[np.argmax(ys)]/10:
                offset = xs[np.argmax(ys)]
            else:
                offset = xs[np.argmin(abs(xs))]
        else:
            offset = xs[np.argmin(abs(xs))]
    
    #Having selected our offset peak, determine its 80%-max width to construct a gaussian which will give the width of our final ppm filter
    peak_widths = sp.signal.peak_widths(ppm_dist[1], peaks, rel_height = .8)
    #This line returns the rounded value of the 80%-max width found by matching the offset to its position in xs, and feeding that index position into peaks - a list of indices, returning the peaks-list index of the xs index. The peaks index corresponds to the peak-widths index, returning the width.
    offset_peak_width = np.round(np.asarray(peak_widths)[0, list(peaks).index(list(peaks)[list(xs).index(offset)])])
    return (offset, offset_peak_width)

def find_rt_duplicates(sum_df):
    proteins=[]
    for name in set(sum_df['name'].values):
        proteins.append(sum_df[sum_df['name']==name])

    count = 0
    duplicate_charges = []
    for protein in proteins:
        if len(protein['charge'].values) != len(set(protein['charge'].values)):
            duplicate_charges.append(protein)

    charge_info = []
    for protein in duplicate_charges:
        name_buffer = []
        for charge in set(protein['charge'].values):
            subdf = protein[protein['charge']==charge]
            if len(subdf)>1:
                dts, rts = subdf["im_mono"].values, subdf["RT"].values
                name_buffer.append([(i, j) for i in dts for j in rts])
        charge_info.append(name_buffer)

    protein_names = dict.fromkeys(sum_df['name'].values)
    protein_rts = dict.fromkeys(sum_df['name'].values)
    rt_matches = dict.fromkeys(sum_df['name'].values)
    for key in protein_names.keys():
        protein_names[key] = sum_df.loc[sum_df['name'] == key]
        protein_rts[key] = protein_names[key][['name',"RT"]].values
        rt_matches[key] = dict.fromkeys(protein_rts[key][:,0])
        for tup in protein_rts[key]:
            rt_cluster = np.array([x[0] for x in protein_rts[key] if x[1] == abs(x[1]-tup[1]) <= 0.2])
            lo_line = [x[0] for x in rt_cluster if x[1] == min(rt_cluster[:, 1])]
            hi_line = [x[0] for x in rt_cluster if x[1] == max(rt_cluster[:, 1])]
            rt_matches[key][tup[0]] = (lo_line + hi_line)

    # Clustering working when hits returns empty
    hits = []
    for key in protein_rts.keys():
        for tup in protein_rts[key]:
            rt_cluster = np.array([x[0] for x in protein_rts[key] if x[1] == abs(x[1]-tup[1]) <= 0.2])
            lo_line = [x[0] for x in rt_cluster if x[1] == min(rt_cluster[:, 1])]
            hi_line = [x[0] for x in rt_cluster if x[1] == max(rt_cluster[:, 1])]
        if len(rt_cluster) > 0:
            hits.append(key)
    return (len(hits)==0, hits)

def apply_cluster_weights(dataframe, dt_weight, rt_weight, mz_weight):
    #applies scoring weights to clustering columns in-place, doesn't return
    dataframe['cluster_im'] = dataframe['im_mono'] / dt_weight
    dataframe['cluster_RT'] = dataframe['RT'] / rt_weight
    dataframe['cluster_mz'] = dataframe['mz_mono'] / mz_weight

### SCRIPT ###

#read IMTBX output file
with open(snakemake.input[0]) as file:
    lines = [x.strip() for x in file.readlines()]
out = []
for i, line in enumerate(lines):
    if line == 'SCAN_START':
        RT = float(lines[i+1].split()[2][3:-1])
        j = i+3
        while lines[j] != 'SCAN_END':
            out.append([float(x) for x in lines[j].split()] + [RT])
            j += 1
            
df = pd.DataFrame(out)
df.columns="mz_mono im_mono ab_mono_peak ab_mono_total mz_top im_top ab_top_peak ab_top_total cluster_peak_count idx_top charge mz_cluster_avg ab_cluster_peak ab_cluster_total cluster_corr noise RT".split()

#buffer df
testq = copy.deepcopy(df)

#read list of all proteins in sample
allseq = pd.read_csv(snakemake.input[1])
allseq['mix'] = [2 for x in range(len(allseq))]
allseq['MW'] = [ProtParam.ProteinAnalysis(seq, monoisotopic=True).molecular_weight() for seq in allseq['sequence']]
allseq['len'] = [len(seq) for seq in allseq['sequence']]

#cluster IMTBX lines corresponding to designed sequence estimates, hardcode values are heuristic weights for clustering, all weights are inverse
apply_cluster_weights(testq, dt_weight = 5.0, rt_weight = 0.6, mz_weight = 0.006)

#create dbscan object, fit, and apply cluster ids to testq lines
db = DBSCAN()
db.fit(testq[['cluster_im','cluster_RT','cluster_mz', 'charge']])
clusters = db.fit_predict(testq[['cluster_im','cluster_RT','cluster_mz', 'charge']])
testq['cluster'] = clusters

#for z in range(3): For visualizing cluster characteristics
#    plotcluster(z)

#average clusters within a ppm window of their suspected proteins
sum_df = cluster_df(testq, ppm = 50)

#generate plot of KDE before ppm correction
kde_plot(sum_df, snakemake.output[1])

#identify major peak of abs_ppm_error clusters, apply correction to all monoisotopic mz values
offset, offset_peak_width = find_offset(sum_df)
if offset > 0:
    testq['mz_mono_fix'] = [x * (1000000 - offset) / (1000000) for x in df['mz_mono']]
    testq['mz_mono_fix_round'] = np.round(testq['mz_mono_fix'].values, 3)
else: 
    testq['mz_mono_fix'] = [x * (1000000 + offset) / (1000000) for x in df['mz_mono']]
    testq['mz_mono_fix_round'] = np.round(testq['mz_mono_fix'].values, 3)

#re-cluster on the adjusted MZ, same weights
apply_cluster_weights(testq, dt_weight = 5.0, rt_weight = 0.6, mz_weight = 0.006)

db = DBSCAN()
db.fit(testq[['cluster_im','cluster_RT','cluster_mz', 'charge']])
clusters = db.fit_predict(testq[['cluster_im','cluster_RT','cluster_mz', 'charge']])
testq['cluster'] = clusters

#re-average clusters to single lines, check for duplicate RTs, save sum_df to outfile
sum_df = cluster_df(testq, ppm=math.ceil(offset_peak_width/2), adjusted=True)

#plot adjusted_mz KDE
kde_plot(sum_df, snakemake.output[2])

#check for duplicate RT-groups THIS MAY BE USELESS TODO
no_duplicates, hits = find_rt_duplicates(sum_df)
print("No rt Duplicates: "+str(no_duplicates))
if not no_duplicates:
    print("DUPLICATES: "+hits)

#send sum_df to snakemake output
sum_df.to_csv(snakemake.output[0], index=False)