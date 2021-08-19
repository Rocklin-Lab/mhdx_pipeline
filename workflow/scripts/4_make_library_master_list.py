"""Example Google style docstrings.

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html

"""
import os
import sys  
import copy
import glob
import yaml
import pymzml
import argparse
import peakutils
import statistics
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from fastdtw import fastdtw
from collections import OrderedDict
from scipy.signal import find_peaks
from matplotlib import pyplot as plt
from scipy.spatial.distance import euclidean
from hdx_limit.core.io import limit_read

matplotlib.use("Agg")


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

    # Finds any rt-group index intersections and resolves.
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
    """Resolves remianing intersections of subset-filtered rt-clusters.

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
    """Searches .mzML for LC-dimension extrema and magnitude.

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
    """Normalizes tic magnitude to 1.

    Args:
        tic (np_array): Chromatogram of Total Ionic Current of an LC-MS run

    Returns:
        tic (np_array): normalized tic

    """
    tic = tic / (np.sum(tic, axis=0) + 1)
    return tic


def gen_stretched_times(tic_file_list, stretched_times_plot_outpath=None):
    """Generate all warp-paths between the undeuterated reference .tic and all others .tic files.

    Args:
        tic_file_list (list of strings): list of paths/to/file.tics where 0th index is reference .tic

    Returns:
        stretched_ts1_times (nested list): all rt-labels stretching the unduterated to later timepoints
        stretched_ts2_times (nested list): all rt-labels stretching later timepoints to the undeuterated
    
    """
    ref_tic_dict = limit_read(tic_file_list[0])
    ref_tic_base_sum = ref_tic_dict['tics_base_sums']
    ref_tic_cumulative_sum = ref_tic_dict['tic_cumulative_sum']
    
    stretched_ts1_times = []
    stretched_ts2_times = []

    # Makes stretched times plot.
    fig, ax = plt.subplots()

    for index, tic_file in enumerate(tic_file_list):

        tic_dict = limit_read(tic_file)
        tic_cumulative_sum = tic_dict['tic_cumulative_sum']
        tic_base_sum = tic_dict['tics_base_sums']

        tic_cumsum = ((tic_cumulative_sum.T / (tic_base_sum + 1)) * ref_tic_base_sum).T

        dist, path = gen_warp_path_for_timepoints(ref_tic_cumulative_sum, tic_cumsum)

        stretched_ts1 = path_to_stretch_times(path, 0)
        stretched_ts2 = path_to_stretch_times(path, 1)

        stretched_ts1_times.append(stretched_ts1)
        stretched_ts2_times.append(stretched_ts2)
        if stretched_times_plot_outpath is not None:
            ax.plot(stretched_ts1, label="tic_file_" + str(index))
            ax.set_ylabel("stretched_ts1_times")
            ax.set_xlabel("index")

    if stretched_times_plot_outpath is not None:
        plt.figure(fig.number)
        plt.legend()
        plt.savefig(stretched_times_plot_outpath)

    return stretched_ts1_times, stretched_ts2_times


def main(names_and_seqs_path,
         undeut_mzml,
         intermediates,
         tics,
         mzml_sum_paths,
         timepoints,
         return_flag=None,
         out_path=None,
         rt_group_cutoff=0.2,
         stretched_times_plot_outpath=None,
         normalization_factors_outpath=None,
         normalization_factors_plot_outpath=None):
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
        out_dict (dict): containing the below key/value pairs.
            library_info (dict): library_info as pandas-ready dict.
            normalization_factors(dict): dict of mzml(i)/mzml[0] intensity ratios with identifying information. 
    
    """
    name_and_seq = pd.read_csv(names_and_seqs_path)

    # If plot is none, function runs without plotting.
    stretched_ts1_times, stretched_ts2_times = gen_stretched_times(tics, stretched_times_plot_outpath=stretched_times_plot_outpath)

    lo_time, hi_time, n_lc_timepoints = set_global_scan_bounds(undeut_mzml)

    # Lists undeuterated dataframes to warp and concatenate.
    undfs = []
    for file in intermediates:
        undfs.append(pd.read_csv(file))

    # Applies warp from provided undeut_mzml to each other undeut mzml, including itself.
    for i in range(len(undfs)):
        undfs[i]["pred_RT"] = [
            pred_time(rt, stretched_ts1_times[i], lo_time, hi_time,
                      n_lc_timepoints) for rt in undfs[i]["RT"]
        ]
        undfs[i]["UN"] = [i for line in undfs[i]["RT"]
                         ]  #apply source index to each line

    # Combines undfs and sort.
    catdf = pd.concat(undfs)
    catdf = catdf.sort_values(["name", "charge", "pred_RT", "abs_ppm"])
    catdf.index = range(len(catdf))

    # Clears duplicate charges close in RT.
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

    # Adds sequences to dataframe.
    catdf["sequence"] = [
        list(name_and_seq.loc[name_and_seq["name"] == catdf.iloc[i]["name"]]
             ["sequence"].values)[0] for i in range(len(catdf))
    ]

    # Applies index after sorting and removing duplicates.
    catdf["idx"] = [i for i in range(len(catdf))]

    # Clusters RT values and renames.
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

    # Drops duplicate charge states, keeps lower error duplicates.
    catdf = catdf.sort_values(["name", "charge", "abs_ppm"])
    catdf = catdf.drop_duplicates(subset=["name", "charge"])

    # Makes rt-group averages weighted by total intensity.
    weighted_avgs = {}
    for name in set(catdf["name"].values):
        weighted_avgs[name] = np.average(
            catdf.loc[catdf["name"]==name]["pred_RT"].values,
            weights=catdf.loc[catdf["name"]==name]["ab_cluster_total"])
    
    # Applies weighted avg to all rt-group members.
    catdf["weighted_average_rt"] = [weighted_avgs[x] for x in catdf["name"].values]

    catdf = catdf.sort_values(["weighted_average_rt", "charge"])
    catdf.index = range(len(catdf))

    # Creates RT_n_m names, where n is the index of the timepoint the source tic came from, and m is the filename index of the tic sourcefile in config[timepoint].
    rt_columns = []
    for i in range(len(timepoints["timepoints"])):
        base = "RT_%s" % i
        if len(timepoints[timepoints["timepoints"][i]]) > 1:
            for j in range(len(timepoints[timepoints["timepoints"][i]])):
                rt_columns.append(base + "_%s" % j)
        else:
            rt_columns.append(base + "_0")

    # Applies warp from provided undeut_mzml RT to each later timepoint RT for each charged species identified.
    for i, stretched in enumerate(stretched_ts2_times):
        catdf[rt_columns[i]] = [
            pred_time(x, stretched, lo_time, hi_time, n_lc_timepoints)
            for x in catdf["pred_RT"]
        ]

    # Determines rt-group average pred-RT-n times from above.
    prev_name = None
    all_tp_mean_preds = [[] for i in range(len(rt_columns))]
    catdf = catdf.sort_values(["weighted_average_rt", "name"])
    for i in range(len(catdf)):
        if catdf.iloc[i]["name"] != prev_name:
            # Get sub frame of rt-group.
            protein_name = catdf.iloc[i]["name"]
            subdf = catdf.loc[catdf["name"] == protein_name]
            # Takes weighted-avg of rt-tp-predictions for all charges in rt-group, if single species group, use species pred-rts as 'mean' stand-ins.
            if len(subdf) > 1:
                name_rt_preds = [
                    np.average(subdf.iloc[:, j].values, weights=catdf.loc[catdf["name"]==protein_name]["ab_cluster_total"])
                    for j in np.arange(-len(rt_columns), 0, 1)
                ]
            else:
                name_rt_preds = subdf.iloc[0, -len(rt_columns):].values
            # Sets avg rt preds for all lines in rt-group
            [[
                all_tp_mean_preds[i].append(name_rt_preds[i])
                for i in range(len(all_tp_mean_preds))
            ]
             for j in range(len(subdf))]

            prev_name = catdf.iloc[i]["name"]
        else:
            pass
    # Sets new columns to give all lines their rt-group RT_n consensus rt-positions.
    for i in range(len(all_tp_mean_preds)):
        catdf["rt_group_mean_" + rt_columns[i]] = all_tp_mean_preds[i]

    ref_mzml_path = [mzml_path for mzml_path in mzml_sum_paths if timepoints[0][0] in mzml_path][0] # Default first undeuterated replicate. 
    ref_sum = float(open(ref_mzml_path, 'r').read())
    # Initializes normalization_factors dict with reference mzml.
    normalization_factors = {"mzml": ["_".join(ref_mzml_path.split("/")[-1].split("_")[:-1])], "sum": [ref_sum], "normalization_factor": [1]}
    for mzml_sum_path in mzml_sum_paths[1:]:
        my_sum = float(open(mzml_sum_path, 'r').read())
        my_mzml = "_".join(mzml_sum_path.split("/")[-1].split("_")[:-1]) #expects path/to/<mzml>_sum.txt
        normalization_factors["mzml"].append(my_mzml)
        normalization_factors["sum"].append(my_sum)
        normalization_factors["normalization_factor"].append(my_sum/ref_sum)

    # Handles output options:
    if out_path is not None:
        catdf.to_json(out_path)

    if normalization_factors_outpath is not None:
        pd.DataFrame.from_dict(normalization_factors).to_csv(normalization_factors_outpath, index=False)

    if normalization_factors_plot_outpath is not None:
        fig = plt.figure()
        fig.suptitle("Normalization Factor Magnitudes")
        ax1 = fig.add_subplot(111)
        ax1.bar(range(len(normalization_factors["mzml"])), normalization_factors["normalization_factor"]) 
        ax1.set(xlabel="TIC .mzML Source", ylabel = "Normalization Factor Magnitude")
        plt.savefig(normalization_factors_plot_outpath)
    
    if return_flag is not None:
        return {"library_info": catdf.to_dict(), "normalization_factors": normalization_factors}


if __name__ == "__main__":

    if "snakemake" in globals():
        names_and_seqs_path = snakemake.input[0]
        open_timepoints = yaml.load(open(snakemake.input[1], "rt"),Loader=yaml.FullLoader)
        undeut_mzml = [fn for fn in snakemake.input if fn.endswith('.mzML.gz')][0]
        tics = [fn for fn in snakemake.input if '.tic' in fn]
        intermediates = [fn for fn in snakemake.input if '_intermediate.csv' in fn]
        mzml_sum_paths = [fn for fn in snakemake.input if '_sum.txt' in fn]
        out_path = snakemake.output[0]
        stretched_times_plot_outpath = snakemake.output[1]
        normalization_factors_outpath = snakemake.output[2]
        normalization_factors_plot_outpath = snakemake.output[3]

        main(names_and_seqs_path=names_and_seqs_path,
             out_path=out_path,
             undeut_mzml=undeut_mzml,
             intermediates=intermediates,
             tics=tics,
             mzml_sum_paths=mzml_sum_paths,
             timepoints=open_timepoints,
             stretched_times_plot_outpath=stretched_times_plot_outpath,
             normalization_factors_outpath=normalization_factors_outpath,
             normalization_factors_plot_outpath=normalization_factors_plot_outpath)

    else:
        parser = argparse.ArgumentParser(
            description=
            "Creates a list of library proteins observed in HDX-LC-IM-MS from imtbx .peaks.isotopes, undeuterated .mzML, and .ims.mz.tic files."
        )
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
            "-v",
            "--mzml_sum_paths",
            nargs="+",
            help="used in snakemake, list of all mzml_sum.txt file paths")
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
        parser.add_argument("-p",
                            "--stretched_times_plot_outpath",
                            help="path/to/stretched_times_plot.png")
        parser.add_argument("-o",
                            "--out_path",
                            help="path/to/library_info.json main output file")
        parser.add_argument("-f",
                            "--normalization_factors_outpath",
                            help="path/to/normalization_factors.csv")
        parser.add_argument("-l",
                            "--normalization_factors_plot_outpath",
                            help="path/to/normalization_factors_plot.png")
        args = parser.parse_args()

        # Generates explicit filenames and open timepoints .yaml.
        if args.mzml_dir is not None and args.undeut_match_string is not None and args.undeut_mzMLs is None:
            args.undeut_mzml = list(
                glob.glob(args.mzml_dir + "*" + args.undeut_match_string + "*" + ".mzML"))
        if args.intermediates_dir is not None and args.intermediates is None:
            args.intermediates = list(
                glob.glob(args.intermediates_dir + "*intermediate.csv"))
        if args.tics_dir is not None and args.tics is None:
            args.tics = list(glob.glob(args.tics_dir + "*.ims.mz.tic.cpickle.zlib"))
        open_timepoints = yaml.load(open(args.timepoints, "rt"),
                                    Loader=yaml.FullLoader)

        main(args.names_and_seqs_path,
             out_path=args.out_path,
             undeut_mzml=args.undeut_mzml,
             intermediates=args.intermediates,
             tics=args.tics,
             mzml_sum_paths=args.mzml_sum_paths,
             timepoints=open_timepoints,
             rt_group_cutoff=args.rt_group_cutoff,
             stretched_times_plot_outpath=args.stretched_times_plot_outpath,
             normalization_factors_outpath=args.normalization_factors_outpath,
             normalization_factors_plot_outpath=args.normalization_factors_plot_outpath)
