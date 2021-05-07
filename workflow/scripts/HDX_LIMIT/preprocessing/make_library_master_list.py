import os
import sys
import copy
import glob
import yaml
import pymzml
import argparse
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

###Definitions###


def path_to_stretch_times(path, to_stretch=0):
    """Applies timewarp to one of the two tics represented in path, 0 to warp the undeuterated-reference to the target, and 1 for the reverse. 

    Args:
        path (list): Dynamic Time Warping least-cost path between two chromatograms
        to_stretch (int): Indicates which path is stretched to the other, in practice 0 is the undeuterated and 1 is a target timepoint

    Returns:
        out_times (list): List of stretched rt-times mapping one tic to the other

    """
    # Applies transformation defined by minimum-cost path from fastdtw to timeseries data
    alt = 1 if to_stretch == 0 else 0

    path = np.array(path)
    out_times = []
    # take the average of the test-equivalent values corresponding to the value at the reference pattern
    for i in range(max(path[:, to_stretch])):
        out_times.append(np.average(path[:, to_stretch][path[:, alt] == i]))
    return np.array(out_times)


def pred_time(rt, stretched_times, lo_time, hi_time, n_lc_timepoints):
    """Converts stretched LC bins to absolute LC retention-time.

    Args:
        rt (float): a point in LC retention-time
        stretched_times (list): Remapped bin labels from path_to_stretch_times
        lo_time (float): lowest scan-time of reference chromatogram
        hi_time (float): highest scan-time of reference chromatogram
        n_lc_timepoints (int): number of LC bins in reference chromatogram

    Returns:
        (float): warped point in LC-RT

    """
    time = int(((rt - lo_time) / (hi_time - lo_time)) * n_lc_timepoints)
    return ((stretched_times[time] / n_lc_timepoints) *
            (hi_time - lo_time)) + lo_time


def rt_cluster(df, name_dict, key, rt_group_cutoff):
    """Groups and filters identified charged-species by rt-distance.

    Parameters:
        df (Pandas DataFrame): df containing all identified charged species
        name_dict (dict): dictionary to be filled with rt-group-member library_info indices 
        key (string): a single library_protein.pdb string
    
    Returns:
        None

    """
    n_df = df.loc[df["name"] == key]
    clusters = [[
        j for j in n_df["idx"].values if
        abs(n_df.loc[n_df["idx"] == i]["pred_RT"].values[0] -
            n_df.loc[n_df["idx"] == j]["pred_RT"].values[0]) < rt_group_cutoff
    ] for i in n_df["idx"].values]
    no_dups = []
    [no_dups.append(lst) for lst in clusters if lst not in no_dups]
    name_dict[key] = subset_filter(no_dups, n_df)


def subset_filter(clusters, n_df):
    """Determines if any rt=cluster is a subset of any other cluster, and removes them. 

    Args:
        clusters (list of lists of ints): List of all clusters
        n_df (Pandas DataFrame): DF of all charged species identified as one library protein 

    Returns:
        final (list of list of ints): mutated input list with all subset rt-groups removed
    
    """
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

    #find any rt-group index intersections and resolve
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
    """Resolve remianing intersections of subset-filtered rt-clusters.

    Args:
        final (list of list  of ints): output of subset_filter, list of lists of all rt-cluster indices
        intersections (list of ints): list of indices in more than one rt-group
        n_df (Pandas DataFrame): DF of all charged species identified as one library protein

    Returns:
        final_copy (list of list of ints): remapped rt-groups with no intersections
    
    """
    final_copy = copy.deepcopy(final)
    [[final_copy[i].discard(j)
      for j in intersections]
     for i in range(len(final_copy))
    ]  # modifies final_copy in place, removing intersection values from each cluster in final_copy
    means = [
        np.mean(n_df.loc[n_df["idx"].isin(list(st))]["pred_RT"].values)
        for st in final_copy
    ]  # generates mean pred_RT of each mutualy exclusive cluster, order maps to final_copy
    dists = [[
        abs(n_df.loc[n_df["idx"] == i]["pred_RT"].values - mean)
        for mean in means
    ]
             for i in intersections
            ]  # outer order maps to intersections, inner maps to final_copy
    [
        final_copy[dists[i].index(min(dists[i]))].add(intersections[i])
        for i in range(len(intersections))
    ]  # adds value from intersection to best-fitting cluster
    return final_copy


def set_global_scan_bounds(mzml):
    """Search .mzML for LC-dimension extrema and magnitude.

    Args:
        mzml (string): path/to/undeuterated.mzML

    Returns:
        lo_time (float): lowest scan-time of reference chromatogram
        hi_time (float): highest scan-time of reference chromatogram
        n_lc_timepoints (int): number of LC bins in reference chromatogram
    
    """
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

    n_lc_timepoints = lc_times - lo_lc_tp

    return np.round(lo_time, 3), np.round(hi_time, 3), n_lc_timepoints


def gen_warp_path_for_timepoints(reference_tic, target_tic):
    """Applies the fast dynamic time-warping algorithm to two provided tics, returns a minimum-cost path to use as a stretching function.

    Args:
        reference_tic (np_array): Undeuterated .tic to be used as reference in warping
        target_tic (np_array): some other .tic to warp to the reference

    Returns:
        distance (float): length of warped path
        path (list): a mapping between the indices of the two tics
    
    """
    distance, path = fastdtw(reference_tic.T,
                             target_tic.T,
                             dist=euclidean,
                             radius=20)
    return distance, path


def norm_tic(tic):
    """Normalize tic magnitude to 1.

    Args:
        tic (np_array): Chromatogram of Total Ionic Current of an LC-MS run

    Returns:
        tic (np_array): normalized tic

    """
    tic = tic / (np.sum(tic, axis=0) + 1)
    return tic


def gen_stretched_times(tic_file_list, plot_path=None):
    """Generate all warp-paths between the undeuterated reference .tic and all others .tic files.

    Args:
        tic_file_list (list of strings): list of paths/to/file.tics where 0th index is reference .tic

    Returns:
        stretched_ts1_times (nested list): all rt-labels stretching the unduterated to later timepoints
        stretched_ts2_times (nested list): all rt-labels stretching later timepoints to the undeuterated
    
    """
    ref_tic = np.loadtxt(tic_file_list[0])
    ref_tic_norm = norm_tic(ref_tic)

    stretched_ts1_times = []
    stretched_ts2_times = []

    # gen a plot
    fig, ax = plt.subplots()

    for index, tic_file in enumerate(tic_file_list):
        tic = np.loadtxt(tic_file)
        tic_norm = norm_tic(tic)

        dist, path = gen_warp_path_for_timepoints(ref_tic_norm, tic_norm)

        stretched_ts1 = path_to_stretch_times(path, 0)
        stretched_ts2 = path_to_stretch_times(path, 1)

        stretched_ts1_times.append(stretched_ts1)
        stretched_ts2_times.append(stretched_ts2)
        if plot_path is not None:
            ax.plot(stretched_ts1, label="tic_file_" + str(index))
            ax.set_ylabel("stretched_ts1_times")
            ax.set_xlabel("index")

    if plot_path is not None:
        plt.legend()
        # save the plot
        plt.savefig(plot_path)

    return stretched_ts1_times, stretched_ts2_times


##########################################################
####################    Operation    s####################
##########################################################


def main(names_and_seqs_path,
         undeut_mzml,
         intermediates,
         tics,
         timepoints,
         return_flag=None,
         out_path=None,
         rt_group_cutoff=0.2,
         plot=None):
    """Generates the master list of library_proteins identified in MS data: library_info.csv.

    Args:
        names_and_seqs_path (string): path/to/names_and_seqs.csv
        undeut_mzml (string): path/to/undeuterated.mzML
        intermediates (list of strings): list of paths to imtbx intermediate files
        tics (list of strings): list of paths to all .tic files
        timepoints (dict): dictionary with 'timepoints' key containing list of hdx timepoints in integer seconds, which are keys mapping to lists of each timepoint's replicate .mzML filenames 
        return_flag (any non-None type): option to return main output in python, for notebook context
        out_path (string): path/to/file for main output library_info.csv
        rt_group_cutoff (float): radius in LC-RT to consider signals a part of an rt-cluster
        plot (any non-None type): path/to/file for stretched time plots

    Returns:
        library_info (dict): Outputs library_info as dict
    
    """
    name_and_seq = pd.read_csv(names_and_seqs_path)

    #if plot is none, function runs without plotting
    stretched_ts1_times, stretched_ts2_times = gen_stretched_times(tics, plot)

    lo_time, hi_time, n_lc_timepoints = set_global_scan_bounds(undeut_mzml)

    # merge UNs - no, appends open dfs to list, why?
    undfs = []
    for file in intermediates:
        undfs.append(pd.read_csv(file))

    # applies warp from provided undeut_mzml to each other undeut mzml, including itself
    for i in range(len(undfs)):
        undfs[i]["pred_RT"] = [
            pred_time(rt, stretched_ts1_times[i], lo_time, hi_time,
                      n_lc_timepoints) for rt in undfs[i]["RT"]
        ]
        undfs[i]["UN"] = [i for line in undfs[i]["RT"]
                         ]  #apply source index to each line

    # combine undfs and sort
    catdf = pd.concat(undfs)
    catdf = catdf.sort_values(["name", "charge", "pred_RT"])
    catdf.index = range(len(catdf))

    # clear duplicate lines TODO: this can be done with DataFrame.drop_duplicates()
    dups = [False]
    for i in range(1, len(catdf)):
        if ((catdf["name"].values[i] == catdf["name"].values[i - 1]) and
            (catdf["charge"].values[i] == catdf["charge"].values[i - 1]) and
            (abs(catdf["pred_RT"].values[i] - catdf["pred_RT"].values[i - 1]) <
             rt_group_cutoff)):
            dups.append(True)
        else:
            dups.append(False)
    catdf["dup"] = dups
    catdf = catdf.query("dup == False")
    catdf["sequence"] = [
        list(name_and_seq.loc[name_and_seq["name"] == catdf.iloc[i]["name"]]
             ["sequence"].values)[0] for i in range(len(catdf))
    ]  # adds sequences to output
    catdf["idx"] = [i for i in range(len(catdf))]

    # cluster RT values and rename
    name_dict = OrderedDict.fromkeys(catdf["name"].values)
    [
        rt_cluster(catdf, name_dict, key, rt_group_cutoff)
        for key in name_dict.keys()
    ]  # TODO possibly automate rt_group cutoff determination in the future

    for key in name_dict.keys():
        for cluster in name_dict[key]:
            mean = np.mean(catdf.iloc[list(cluster)]["pred_RT"].values)
            for line in list(cluster):
                catdf.iat[line, 0] = catdf.iloc[line]["name"] + "_" + str(
                    round(mean, 5))

    # Make rt-group averages weighted by total intensity.
    weighted_avgs = {}
    for name in set(catdf["name"].values):
        weighted_avgs[name] = np.average(
            catdf.loc[catdf["name"]==name]["pred_RT"].values,
            weights=catdf.loc[catdf["name"]==name]["ab_cluster_total"])
    #apply weighted avg to all rt-group members
    catdf["weighted_average_rt"] = [weighted_avgs[x] for x in catdf["name"].values]

    catdf = catdf.sort_values(["weighted_average_rt", "charge"])
    catdf.index = range(len(catdf))

    # Create RT_n_m names, where n is the index of the timepoint the source tic came from, and m is the filename index of the tic sourcefile in config[timepoint]
    rt_columns = []
    for i in range(len(timepoints["timepoints"])):
        base = "RT_%s" % i
        if len(timepoints[timepoints["timepoints"][i]]) > 1:
            for j in range(len(timepoints[timepoints["timepoints"][i]])):
                rt_columns.append(base + "_%s" % j)
        else:
            rt_columns.append(base + "_0")

    # apply warp from provided undeut_mzml RT to each later timepoint RT for each charged species identified
    for i, stretched in enumerate(stretched_ts2_times):
        catdf[rt_columns[i]] = [
            pred_time(x, stretched, lo_time, hi_time, n_lc_timepoints)
            for x in catdf["pred_RT"]
        ]

    # determine rt-group average pred-RT-n times from above
    prev_name = None
    all_tp_mean_preds = [[] for i in range(len(rt_columns))]
    catdf = catdf.sort_values(["weighted_average_rt", "name"])
    for i in range(len(catdf)):
        if catdf.iloc[i]["name"] != prev_name:
            # get sub frame of rt-group
            protein_name = catdf.iloc[i]["name"]
            subdf = catdf.loc[catdf["name"] == protein_name]
            # take weighted-avg of rt-tp-predictions for all charges in rt-group, if single species group, use species pred-rts as 'mean' stand-ins
            if len(subdf) > 1:
                name_rt_preds = [
                    np.average(subdf.iloc[:, j].values, weights=catdf.loc[catdf["name"]==protein_name]["ab_cluster_total"])
                    for j in np.arange(-len(rt_columns), 0, 1)
                ]
            else:
                name_rt_preds = subdf.iloc[0, -len(rt_columns):].values
            # Set avg rt preds for all lines in rt-group
            [[
                all_tp_mean_preds[i].append(name_rt_preds[i])
                for i in range(len(all_tp_mean_preds))
            ]
             for j in range(len(subdf))]
            # set prev_name to current name
            prev_name = catdf.iloc[i]["name"]
        else:
            pass
    # set new columns to give all lines their rt-group RT_n consensus rt-positions
    for i in range(len(all_tp_mean_preds)):
        catdf["rt_group_mean_" + rt_columns[i]] = all_tp_mean_preds[i]

    if out_path is not None:
        catdf.to_csv(out_path)

    if return_flag is not None:
        return catdf.to_dict()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=
        "Creates a list of library proteins observed in HDX-LC-IM-MS from imtbx .peaks.isotopes, undeuterated .mzML, and .ims.mz.tic files."
    )
    # inputs
    parser.add_argument(
        "names_and_seqs_path",
        help="path/to/file .csv of protein library names and sequences")
    parser.add_argument("-m",
                        "--mzml_dir",
                        help="path/to/dir/ containing undeuterated .mzML files")
    parser.add_argument(
        "-s",
        "--undeut_match_string",
        help="unique part of undeuterated mzML filename to be used in matching")
    parser.add_argument("-i",
                        "--intermediates_dir",
                        help="path/to/dir/ containing intermediate imtbx files")
    parser.add_argument("-t",
                        "--tics_dir",
                        help="path/to/dir/ containing .ims.mz.tic files")
    parser.add_argument("-n",
                        "--undeut_mzml",
                        help="path/to/file, one undeuterated .mzML")
    parser.add_argument(
        "-j",
        "--intermediates",
        nargs="+",
        help="used in snakemake, list of all imtbx intermediate file paths")
    parser.add_argument(
        "-u",
        "--tics",
        nargs="+",
        help="used in snakemake, list of all .imx.mz.tic file paths")
    parser.add_argument(
        "-e",
        "--timepoints",
        required=True,
        help=
        "path/to/.yaml file with snakemake.config timepoints and .mzML filenames by timepoint"
    )
    parser.add_argument(
        "-c",
        "--rt_group_cutoff",
        default=0.2,
        type=float,
        help=
        "control value for creation of RT-groups, maximum rt-distance between same-mass isotope clusters"
    )
    # outputs
    parser.add_argument("-p",
                        "--plot",
                        help="path/to/stretched_times_plots.png")
    parser.add_argument("-o",
                        "--out_path",
                        help="path/to/library_info.csv main output file")

    args = parser.parse_args()

    #Generate explicit filenames and open timepoints .yaml
    if args.mzml_dir is not None and args.undeut_match_string is not None and args.undeut_mzMLs is None:
        args.undeut_mzml = list(
            glob.glob(args.mzml_dir + "*" + args.undeut_match_string + "*" +
                      ".mzML"))
    if args.intermediates_dir is not None and args.intermediates is None:
        args.intermediates = list(
            glob.glob(args.intermediates_dir + "*intermediate.csv"))
    if args.tics_dir is not None and args.tics is None:
        args.tics = list(glob.glob(args.tics_dir + "*.ims.mz.tic"))
    open_timepoints = yaml.load(open(args.timepoints, "rt"),
                                Loader=yaml.FullLoader)

    main(args.names_and_seqs_path,
         out_path=args.out_path,
         undeut_mzml=args.undeut_mzml,
         intermediates=args.intermediates,
         tics=args.tics,
         timepoints=open_timepoints,
         rt_group_cutoff=args.rt_group_cutoff,
         plot=args.plot)
