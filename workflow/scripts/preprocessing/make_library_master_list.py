import copy
import pymzml
import peakutils
import statistics
import numpy as np
import pandas as pd
import seaborn as sns
from fastdtw import fastdtw
from collections import OrderedDict
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.spatial.distance import euclidean

import matplotlib

matplotlib.use("Agg")

import ipdb


###Definitions###


def path_to_stretch_times(path, to_stretch=0):
    # Applies transformation defined by minimum-cost path from fastdtw to timeseries data
    alt = 1 if to_stretch == 0 else 0

    path = np.array(path)
    out_times = []
    for i in range(max(path[:, to_stretch])):
        # take the average of the test-equivalent values corresponding to the value at the reference pattern
        out_times.append(np.average(path[:, to_stretch][path[:, alt] == i]))
    return np.array(out_times)


def pred_time(x, stretched_times, lo_time, hi_time, lc_timepoints):
    time = int(
        ((x - lo_time) / (hi_time - lo_time)) * lc_timepoints
    )  ###previously hardcoded as lo=3, hi=17, diff=987 (987 not correct, but not really harmful, ~1-2% error)
    return ((stretched_times[time] / lc_timepoints) * (hi_time - lo_time)) + lo_time


def cluster(df, name_dict, key, RT_cutoff):
    n_df = df.loc[df["name"] == key]
    clusters = [
        [
            j
            for j in n_df["idx"].values
            if abs(
                n_df.loc[n_df["idx"] == i]["pred_RT"].values[0]
                - n_df.loc[n_df["idx"] == j]["pred_RT"].values[0]
            )
            < RT_cutoff
        ]
        for i in n_df["idx"].values
    ]
    no_dups = []
    [no_dups.append(lst) for lst in clusters if lst not in no_dups]
    name_dict[key] = subset_filter(no_dups, n_df)


def subset_filter(clusters, n_df):
    sets = [set(cluster) for cluster in clusters]
    final = []
    for i in range(len(sets)):
        sup = True
        for j in range(len(sets)):
            if i != j:
                if sets[i].issubset(sets[j]):
                    sup = False
        if sup:
            final.append(sets[i])

    intersections = []
    for s1 in final:
        for i in s1:
            for s2 in final:
                if i in s2 and s1 != s2:
                    intersections.append(i)
    intersections = list(set(intersections))

    if len(intersections) > 0:
        return intersection_filter(final, intersections, n_df)
    else:
        return final


def intersection_filter(final, intersections, n_df):
    final_copy = copy.deepcopy(final)
    [
        [final_copy[i].discard(j) for j in intersections]
        for i in range(len(final_copy))
    ]  # modifies final_copy in place, removing intersection values from each cluster in final_copy
    means = [
        np.mean(n_df.loc[n_df["idx"].isin(list(st))]["pred_RT"].values)
        for st in final_copy
    ]  # generates mean pred_RT of each mutualy exclusive cluster, order maps to final_copy
    dists = [
        [abs(n_df.loc[n_df["idx"] == i]["pred_RT"].values - mean) for mean in means]
        for i in intersections
    ]  # outer order maps to intersections, inner maps to final_copy
    [
        final_copy[dists[i].index(min(dists[i]))].add(intersections[i])
        for i in range(len(intersections))
    ]  # adds value from intersection to best-fitting cluster
    return final_copy


def set_global_scan_bounds(mzml):
    run = pymzml.run.Reader(mzml)
    n_scans = run.get_spectrum_count()
    lc_times = int(n_scans / 200)
    last = 0
    no_lo = True
    for spectrum in run:
        if spectrum.index % 200 == 0:
            time = spectrum.scan_time_in_minutes()
            if abs(time - last) > 0.005:
                if no_lo:
                    lo_time = spectrum.scan_time_in_minutes()
                    lo_lc_tp = int(spectrum.index // 200)
                no_lo = False
            if int(spectrum.index // 200) == lc_times - 1:
                hi_time = spectrum.scan_time_in_minutes()
            last = time

    lc_timepoints = lc_times - lo_lc_tp

    return np.round(lo_time, 3), np.round(hi_time, 3), lc_timepoints


def gen_warp_path_for_timepoints(reference_tics, tics):
    distance, path = fastdtw(reference_tics.T, tics.T, dist=euclidean, radius=20)
    return distance, path


def norm_tics(tics):
    tics = tics / (np.sum(tics, axis=0) + 1)
    return tics


def gen_stretched_times(tic_file_list):
    ref_tics = np.loadtxt(tic_file_list[0])
    ref_tics_norm = norm_tics(ref_tics)

    stretched_ts1_times = []
    stretched_ts2_times = []

    # gen a plot
    fig, ax = plt.subplots()

    for index, tic_file in enumerate(tic_file_list):
        tics = np.loadtxt(tic_file)
        tics_norm = norm_tics(tics)

        dist, path = gen_warp_path_for_timepoints(ref_tics_norm, tics_norm)

        stretched_ts1 = path_to_stretch_times(path, 0)
        stretched_ts2 = path_to_stretch_times(path, 1)

        stretched_ts1_times.append(stretched_ts1)
        stretched_ts2_times.append(stretched_ts2)

        ax.plot(stretched_ts1, label="tic_file_" + str(index))
        ax.set_ylabel("stretched_ts1_times")
        ax.set_xlabel("index")

    plt.legend()
    # save the plot
    plt.savefig(snakemake.output[0])

    return stretched_ts1_times, stretched_ts2_times


##########################################################
####################    Operation    s####################
##########################################################

ins = snakemake.input

name_and_seq = pd.read_csv(ins.pop(0))  # name and seq list comes first
mzml = ins.pop(0)
csvs = [fn for fn in ins if ".csv" in fn]

tics_filelist = [fn for fn in ins if ".ims.mz.tic" in fn]

stretched_ts1_times, stretched_ts2_times = gen_stretched_times(tics_filelist)

lo_time, hi_time, lc_timepoints = set_global_scan_bounds(mzml)

# merge UNs
undfs = []
for file in csvs:
    undfs.append(pd.read_csv(file))

# inputs the transformed test-pattern value for each undeuterated imtbx csv as pred_RT
for i in range(len(undfs)):
    undfs[i]["pred_RT"] = [
        pred_time(x, stretched_ts1_times[i], lo_time, hi_time, lc_timepoints)
        for x in undfs[i]["RT"]
    ]
    undfs[i]["UN"] = [i for x in undfs[i]["RT"]]

# combine undfs and sort
catdf = pd.concat(undfs)
catdf = catdf.sort_values(["name", "charge", "pred_RT"])
catdf.index = range(len(catdf))

# clear duplicate lines
dups = [False]
for i in range(1, len(catdf)):
    if (
        (catdf["name"].values[i] == catdf["name"].values[i - 1])
        and (catdf["charge"].values[i] == catdf["charge"].values[i - 1])
        and (
            abs(catdf["pred_RT"].values[i] - catdf["pred_RT"].values[i - 1])
            < snakemake.config["rt_group_cutoff"]
        )
    ):  # ensures no duplicate charge-states make it into same rt-group
        dups.append(True)
    else:
        dups.append(False)
catdf["dup"] = dups
catdf = catdf.query("dup == False")
catdf["sequence"] = [
    list(
        name_and_seq.loc[name_and_seq["name"] == catdf.iloc[i]["name"]][
            "sequence"
        ].values
    )[0]
    for i in range(len(catdf))
]  # adds sequences to output
catdf["idx"] = [i for i in range(len(catdf))]

# cluster RT values and rename
name_dict = OrderedDict.fromkeys(catdf["name"].values)
[
    cluster(catdf, name_dict, key, snakemake.config["rt_group_cutoff"])
    for key in name_dict.keys()
]  # TODO possibly automate rt_group cutoff determination in the future

for key in name_dict.keys():
    for cluster in name_dict[key]:
        mean = np.mean(catdf.iloc[list(cluster)]["pred_RT"].values)
        for line in list(cluster):
            catdf.iat[line, 0] = catdf.iloc[line]["name"] + "_" + str(round(mean, 5))

# stick each 'rt_group' entry with a median RT
med_RTs = {}
for name in set(catdf["name"].values):
    med_RTs[name] = np.median(catdf.query('name == "%s"' % name)["pred_RT"].values)
catdf["med_RT"] = [med_RTs[x] for x in catdf["name"].values]

catdf = catdf.sort_values(["med_RT", "charge"])
catdf.index = range(len(catdf))

# Create RT_n_m names, where n is the index of the timepoint the source tic came from, and m is the filename index of the tic sourcefile in config[timepoint]
rt_columns = []
for i in range(len(snakemake.config["timepoints"])):
    base = "RT_%s" % i
    if len(snakemake.config[snakemake.config["timepoints"][i]]) > 1:
        for j in range(len(snakemake.config[snakemake.config["timepoints"][i]])):
            rt_columns.append(base + "_%s" % j)
    else:
        rt_columns.append(base + "_0")

for i, stretched in enumerate(stretched_ts2_times):
    catdf[rt_columns[i]] = [
        pred_time(x, stretched, lo_time, hi_time, lc_timepoints)
        for x in catdf["pred_RT"]
    ]

# determine rt-group average pred-RT-n times from above
prev_name = None
all_tp_mean_preds = [[] for i in range(len(rt_columns))]
catdf = catdf.sort_values(["med_RT", "name"])
for i in range(len(catdf)):
    if catdf.iloc[i]["name"] != prev_name:
        # get sub frame of rt-group
        subdf = catdf.loc[catdf["name"] == catdf.iloc[i]["name"]]
        # take means of rt-tp-predictions for all charges in rt-group, if single species group, use species pred-rts as 'mean' stand-ins
        if len(subdf) > 1:
            name_rt_preds = [
                np.mean(subdf.iloc[:, i].values)
                for i in np.arange(-len(rt_columns), 0, 1)
            ]
        else:
            name_rt_preds = subdf.iloc[0, -len(rt_columns) :].values
        # Set avg rt preds for all lines in rt-group
        [
            [
                all_tp_mean_preds[i].append(name_rt_preds[i])
                for i in range(len(all_tp_mean_preds))
            ]
            for j in range(len(subdf))
        ]
        # set prev_name to current name
        prev_name = catdf.iloc[i]["name"]
    else:
        pass
# set new columns to give all lines their rt-group RT_n consensus rt-positions
for i in range(len(all_tp_mean_preds)):
    catdf["rt_group_mean_" + rt_columns[i]] = all_tp_mean_preds[i]

catdf.to_csv(snakemake.output[1])
