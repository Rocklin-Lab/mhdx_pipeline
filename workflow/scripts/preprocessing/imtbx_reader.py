import importlib.util

hxtools_spec = importlib.util.spec_from_file_location(
    "hxtools.py", "workflow/scripts/auxiliary/hxtools.py"
)
molmass_spec = importlib.util.spec_from_file_location(
    "molmass.py", "workflow/scripts/auxiliary/molmass.py"
)
hxtools = importlib.util.module_from_spec(hxtools_spec)
molmass = importlib.util.module_from_spec(molmass_spec)
hxtools_spec.loader.exec_module(hxtools)
molmass_spec.loader.exec_module(molmass)

import sys
import time
import copy
import math
import argparse
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
import pickle as pk

matplotlib.use("Agg")

import ipdb


### DEFINITIONS ###


def plotcluster(i=0):
    """ Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value

   """
    plt.figure(figsize=(16, 3))
    plt.subplot(141)
    plt.plot(testq[clusters == i]["RT"], testq[clusters == i]["mz_mono"])
    plt.title(
        "mz range = %.4f"
        % (max(testq[clusters == i]["mz_mono"]) - min(testq[clusters == i]["mz_mono"]))
    )

    plt.subplot(142)
    plt.plot(testq[clusters == i]["RT"], testq[clusters == i]["cluster_im"])
    plt.title(
        "im range = %.1f"
        % (
            max(testq[clusters == i]["cluster_im"])
            - min(testq[clusters == i]["cluster_im"])
        )
    )

    plt.subplot(143)
    plt.plot(testq[clusters == i]["RT"], testq[clusters == i]["ab_cluster_total"])

    plt.subplot(144)
    plt.plot(testq[clusters == i]["RT"], testq[clusters == i]["cluster_corr"])
    plt.savefig("../../plots/" + "_RT_cluster_plots.png")


def kde_plot(sum_df, outpath):
    """ Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value

   """
    mykde = gaussian_kde(sum_df["ppm"])
    xs = np.linspace(-50, 50, 10000)
    xs[np.argmax(mykde.evaluate(xs))]
    sns.distplot(sum_df["ppm"])
    plt.xlim([-50, 50])
    plt.grid()
    plt.plot(xs, mykde(xs))
    plt.savefig(outpath)


# mix doesn't serve any purpose in the pipeline TODO
def getnear(x, charge=None, mix=None, ppm=50):
""" Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value
    """
    subdf = allseq
    if mix != None:
        subdf = allseq[allseq["mix"] == mix]
    if charge != None:
        low, high = (
            ((x * charge) - (1.00727 * charge)) * ((1000000 - ppm) / 1000000),
            ((x * charge) - (1.00727 * charge)) * ((1000000 + ppm) / 1000000),
        )
        mlow, mhigh = allseq["MW"] > low, allseq["MW"] < high
        tempdf = allseq[mlow & mhigh].sort_values("MW")[
            ["MW", "mix", "name", "len", "sequence"]
        ]
        tempdf["plus%s" % int(charge)] = [
            (q + (1.00727 * charge)) / charge for q in tempdf["MW"]
        ]
        tempdf["ppm"] = [
            "%.1f" % ((1.0 - (q / x)) * 1000000) for q in tempdf["plus%s" % int(charge)]
        ]
        tempdf["abs_ppm"] = [
            np.abs(((1.0 - (q / x)) * 1000000)) for q in tempdf["plus%s" % int(charge)]
        ]
        return tempdf[
            [
                "plus%s" % int(charge),
                "ppm",
                "abs_ppm",
                "MW",
                "mix",
                "name",
                "len",
                "sequence",
            ]
        ]
    else:
        low, high = x - window, x + window
        mlow, mhigh = allseq["MW"] > low, allseq["MW"] < high
        tempdf = subdf[mlow & mhigh].sort_values("MW")[
            ["MW", "mix", "name", "len", "sequence"]
        ]
        return tempdf


def cluster_df_hq_signals(
    testq, ppm=50, intensity_threshold=1e4, cluster_correlation=0.99, adjusted=False
):
""" Cluster high quality mz signals based on intensities and cluster correlation
    :param testq: dataframe from imtbx
    :param ppm: ppm error to include for mz signals
    :param intensity_threshold: minimum intensity value required for mz signals
    :param cluster_correlation: cluster correlation from imtbx. higher correlation means better isotopic distribution
    :param adjusted: Boolean to indicate if mz signals have already been corrected
    :return: clustered mz signal
    """

    hq_dataframe = testq[
        (testq["cluster_corr"] > cluster_correlation)
        & (testq["ab_cluster_total"] > (intensity_threshold))
    ]

    sum_data = []
    for c in range(0, max(hq_dataframe["cluster"]) + 1):

        cluster_df = hq_dataframe[hq_dataframe["cluster"] == c]

        if (
            len(cluster_df["mz_mono"]) > 1
        ):  # ask Wes why isn't this set to greater than 1?

            charge = np.median(cluster_df["charge"])
            if adjusted:
                mz = np.average(
                    cluster_df["mz_mono_fix_round"],
                    weights=cluster_df["ab_cluster_total"],
                )
            else:
                mz = np.average(
                    cluster_df["mz_mono"], weights=cluster_df["ab_cluster_total"]
                )
            RT = np.average(cluster_df["RT"], weights=cluster_df["ab_cluster_total"])
            im = np.average(
                cluster_df["im_mono"], weights=cluster_df["ab_cluster_total"]
            )

            near = getnear(mz, charge=charge, mix=2, ppm=ppm)

            if len(near) > 0:
                sum_data.append(
                    [
                        near["name"].values[0],
                        RT,
                        im,
                        sum(cluster_df["ab_cluster_total"]),
                        near["MW"].values[0],
                        charge,
                        near["plus%s" % int(charge)].values[0],
                        mz,
                        near["ppm"].values[0],
                        near["abs_ppm"].values[0],
                        c,
                    ]
                )
            if len(near) > 1:
                display(near)
    sum_df = pd.DataFrame(sum_data)
    sum_df.columns = [
        "name",
        "RT",
        "im_mono",
        "ab_cluster_total",
        "MW",
        "charge",
        "expect_mz",
        "obs_mz",
        "ppm",
        "abs_ppm",
        "cluster",
    ]
    sum_df["ppm"] = [float(x) for x in sum_df["ppm"]]
    return sum_df


def calc_mz_ppm_error(obs_mz, thr_mz):
""" Calculate mz ppm error
    :param obs_mz: observed mz values from the experiment
    :param thr_mz: theoreteical or expected mz values based on chemical composition
    :return: ppm error
    """
    ppm_err = 1e6 * (obs_mz - thr_mz) / thr_mz
    return ppm_err


def gen_calib_dict(polyfit_bool=False, **args):
""" Generate calibration dictionary with keywords.
    :param polyfit_bool: True or False
    :param args: arguements
    :return: calibration dictionary
    """

    calib_dict = dict()

    calib_dict["polyfit_bool"] = polyfit_bool

    if polyfit_bool:
        for param, value in args.items():
            calib_dict[param] = value

    return calib_dict


def gen_mz_ppm_error_calib_polyfit(obs_mz, thr_mz, polyfit_deg=1):
""" Use polyfit to generate a function to correlate observed and theoretical mz values. The function is used as calibration
    for the mz values. User can specify the degree of the polyfit
    :param obs_mz: observed mz values
    :param thr_mz: theoretical mz values
    :param polyfit_deg: degree for polynomial fit
    :return: dictionary containing the dataset used for calibration, polyfit, and ppm error before and after calibration
    """

    polyfit_coeffs = np.polyfit(x=obs_mz, y=thr_mz, deg=polyfit_deg)

    obs_mz_corr = apply_polyfit_cal_mz(polyfit_coeffs, obs_mz)

    ppm_error_before_corr = calc_mz_ppm_error(obs_mz, thr_mz)
    ppm_error_after_corr = calc_mz_ppm_error(obs_mz_corr, thr_mz)

    cal_dict = gen_calib_dict(
        polyfit_bool=True,
        thr_mz=thr_mz,
        obs_mz=obs_mz,
        polyfit_coeffs=polyfit_coeffs,
        polyfit_deg=polyfit_deg,
        obs_mz_corr=obs_mz_corr,
        ppm_error_before_corr=ppm_error_before_corr,
        ppm_error_after_corr=ppm_error_after_corr,
    )

    return cal_dict


def apply_polyfit_cal_mz(polyfit_coeffs, mz):
""" Apply polyfit coeff to transform the mz values
    :param polyfit_coeffs: polyfit coefficients
    :param mz: mz values
    :return: transformed mz values
    """
    mz_corr = np.polyval(polyfit_coeffs, mz)
    return mz_corr


def gen_mz_error_calib_output(
    testq,
    calib_pk_fpath,
    polyfit_degree=1,
    ppm_tol=50,
    int_tol=1e4,
    cluster_corr_tol=0.99,
):
""" Generate calibration using the dataframe from imtbx
    :param testq: dataframe from imtbx
    :param calib_pk_fpath: pickle filepath to save calibration information
    :param polyfit_degree: polyfit degree
    :param ppm_tol: ppm tolerance for selecting mz signals for calibration
    :param int_tol: intensity tolerance for selecting mz signals for calibration
    :param cluster_corr_tol: cluster correlation tolerance for selecting mz signals for calibration
    :return: calibration dictionary.
    """

    # generate high quality cluster mz signals
    cluster_hq_df = cluster_df_hq_signals(
        testq=testq,
        ppm=ppm_tol,
        intensity_threshold=int_tol,
        cluster_correlation=cluster_corr_tol,
    )

    # generate calibration dictionary
    calib_dict = gen_mz_ppm_error_calib_polyfit(
        obs_mz=cluster_hq_df["obs_mz"].values,
        thr_mz=cluster_hq_df["expect_mz"].values,
        polyfit_deg=polyfit_degree,
    )

    # save calibration dictionary for further use
    save_pickle_object(calib_dict, calib_pk_fpath)

    return calib_dict


def save_pickle_object(object, fpath):
""" Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value
    """
    with open(fpath, "wb") as outfile:
        pk.dump(object, outfile)


def cluster_df(testq, ppm=50, adjusted=False):
    sum_data = []
    for c in range(0, max(testq["cluster"]) + 1):
        cluster_df = testq[testq["cluster"] == c]
        charge = np.median(cluster_df["charge"])
        if adjusted:
            mz = np.average(
                cluster_df["mz_mono_fix_round"], weights=cluster_df["ab_cluster_total"]
            )
        else:
            mz = np.average(
                cluster_df["mz_mono"], weights=cluster_df["ab_cluster_total"]
            )
        RT = np.average(cluster_df["RT"], weights=cluster_df["ab_cluster_total"])
        im = np.average(cluster_df["im_mono"], weights=cluster_df["ab_cluster_total"])

        near = getnear(mz, charge=charge, mix=2, ppm=ppm)

        if len(near) > 0:
            sum_data.append(
                [
                    near["name"].values[0],
                    RT,
                    im,
                    sum(cluster_df["ab_cluster_total"]),
                    near["MW"].values[0],
                    charge,
                    near["plus%s" % int(charge)].values[0],
                    mz,
                    near["ppm"].values[0],
                    near["abs_ppm"].values[0],
                    c,
                ]
            )
        if len(near) > 1:
            display(near)
    sum_df = pd.DataFrame(sum_data)
    sum_df.columns = [
        "name",
        "RT",
        "im_mono",
        "ab_cluster_total",
        "MW",
        "charge",
        "expect_mz",
        "obs_mz",
        "ppm",
        "abs_ppm",
        "cluster",
    ]
    sum_df["ppm"] = [float(x) for x in sum_df["ppm"]]
    return sum_df


def find_offset(sum_df):
""" Returns suspected systemic ppm error and width of poi peak of run data from sum_df.
    Assumes protein of interest within +/- 50ppm of 0ppm.
    Selects closest peak to 0 if sufficiently prominent.
    """
    # Maybe make this save the distplot too.
    ppm_dist = sns.distplot(sum_df["ppm"].values).get_lines()[0].get_data()
    peaks = sp.signal.find_peaks(ppm_dist[1])[0]
    xs, ys = ppm_dist[0][peaks], ppm_dist[1][peaks]
    # If lowest ppm peak is also highest frequency within window, we hope this will be the common case in our 50 ppm initial window
    try:
        xs[np.argmin(abs(xs))] == ys[np.argmax(ys)]
    except:
        ipdb.set_trace()
    if xs[np.argmin(abs(xs))] == ys[np.argmax(ys)]:
        return xs[np.argmin(abs(xs))]
        # Lowest ppm peak is not most prominent, determine relative height of lowest ppm peak
    else:
        # If peak closer to zero is less than half the height of the more prominent peak, check larger peak's ppm
        if ys[np.argmin(abs(xs))] < ys[np.argmax(ys)] / 2:
            # If most prominent peak is heuristically close to 0, or lowest ppm peak is relatively very small (10% of major peak): use big peak
            if (
                xs[np.argmax(ys)] < 25
                or ys[np.argmin(abs(xs))] < ys[np.argmax(ys)] / 10
            ):
                offset = xs[np.argmax(ys)]
            else:
                offset = xs[np.argmin(abs(xs))]
        else:
            offset = xs[np.argmin(abs(xs))]

    # Having selected our offset peak, determine its 80%-max width to construct a gaussian which will give the width of our final ppm filter
    peak_widths = sp.signal.peak_widths(ppm_dist[1], peaks, rel_height=0.8)
    # This line returns the rounded value of the 80%-max width found by matching the offset to its position in xs, and feeding that index position into peaks - a list of indices, returning the peaks-list index of the xs index. The peaks index corresponds to the peak-widths index, returning the width.
    offset_peak_width = np.round(
        np.asarray(peak_widths)[
            0, list(peaks).index(list(peaks)[list(xs).index(offset)])
        ]
    )
    return (offset, offset_peak_width)


def find_rt_duplicates(sum_df):
""" Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value
    """
    proteins = []
    for name in set(sum_df["name"].values):
        proteins.append(sum_df[sum_df["name"] == name])

    count = 0
    duplicate_charges = []
    for protein in proteins:
        if len(protein["charge"].values) != len(set(protein["charge"].values)):
            duplicate_charges.append(protein)

    charge_info = []
    for protein in duplicate_charges:
        name_buffer = []
        for charge in set(protein["charge"].values):
            subdf = protein[protein["charge"] == charge]
            if len(subdf) > 1:
                dts, rts = subdf["im_mono"].values, subdf["RT"].values
                name_buffer.append([(i, j) for i in dts for j in rts])
        charge_info.append(name_buffer)

    protein_names = dict.fromkeys(sum_df["name"].values)
    protein_rts = dict.fromkeys(sum_df["name"].values)
    rt_matches = dict.fromkeys(sum_df["name"].values)
    for key in protein_names.keys():
        protein_names[key] = sum_df.loc[sum_df["name"] == key]
        protein_rts[key] = protein_names[key][["name", "RT"]].values
        rt_matches[key] = dict.fromkeys(protein_rts[key][:, 0])
        for tup in protein_rts[key]:
            rt_cluster = np.array(
                [x[0] for x in protein_rts[key] if x[1] == abs(x[1] - tup[1]) <= 0.2]
            )
            lo_line = [x[0] for x in rt_cluster if x[1] == min(rt_cluster[:, 1])]
            hi_line = [x[0] for x in rt_cluster if x[1] == max(rt_cluster[:, 1])]
            rt_matches[key][tup[0]] = lo_line + hi_line

    # Clustering working when hits returns empty
    hits = []
    for key in protein_rts.keys():
        for tup in protein_rts[key]:
            rt_cluster = np.array(
                [x[0] for x in protein_rts[key] if x[1] == abs(x[1] - tup[1]) <= 0.2]
            )
            lo_line = [x[0] for x in rt_cluster if x[1] == min(rt_cluster[:, 1])]
            hi_line = [x[0] for x in rt_cluster if x[1] == max(rt_cluster[:, 1])]
        if len(rt_cluster) > 0:
            hits.append(key)
    return (len(hits) == 0, hits)


def apply_cluster_weights(dataframe, dt_weight, rt_weight, mz_weight):
""" Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value
    """
    # TODO This is not great style, this should accept 3 lists and 3 weights and return 3 new lists
    # applies scoring weights to clustering columns in-place, doesn't return
    dataframe["cluster_im"] = dataframe["im_mono"] / dt_weight
    dataframe["cluster_RT"] = dataframe["RT"] / rt_weight
    dataframe["cluster_mz"] = dataframe["mz_mono"] / mz_weight


### SCRIPT ###
def main(isotopes_path, names_and_seqs_path, intermediate_out_path, plot=None, original_mz_kde_path=None, adjusted_mz_kde_path=None, calibration_outpath=None, polyfit_deg=None, ppm_tolerance=None, intensity_tolerance=None, cluster_corr_tolerance=None, ppm_refilter=None):
""" Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value
    """
    
    # read IMTBX output file
    with open(isotopes_path) as file:
        lines = [x.strip() for x in file.readlines()]
    out = []
    for i, line in enumerate(lines):
        if line == "SCAN_START":
            RT = float(lines[i + 1].split()[2][3:-1])
            j = i + 3
            while lines[j] != "SCAN_END":
                out.append([float(x) for x in lines[j].split()] + [RT])
                j += 1

    df = pd.DataFrame(out)
    df.columns = "mz_mono im_mono ab_mono_peak ab_mono_total mz_top im_top ab_top_peak ab_top_total cluster_peak_count idx_top charge mz_cluster_avg ab_cluster_peak ab_cluster_total cluster_corr noise RT".split()

    # buffer df
    testq = copy.deepcopy(df)

    # read list of all proteins in sample
    allseq = pd.read_csv(names_and_seqs_path)
    allseq["mix"] = [2 for x in range(len(allseq))]
    allseq["MW"] = [
        ProtParam.ProteinAnalysis(seq, monoisotopic=True).molecular_weight()
        for seq in allseq["sequence"]
    ]
    allseq["len"] = [len(seq) for seq in allseq["sequence"]]

    # cluster IMTBX lines corresponding to designed sequence estimates, hardcode values are heuristic weights for clustering, all weights are inverse
    apply_cluster_weights(testq, dt_weight=5.0, rt_weight=0.6, mz_weight=0.006)

    # create dbscan object, fit, and apply cluster ids to testq lines
    db = DBSCAN()
    db.fit(testq[["cluster_im", "cluster_RT", "cluster_mz", "charge"]])
    clusters = db.fit_predict(testq[["cluster_im", "cluster_RT", "cluster_mz", "charge"]])
    testq["cluster"] = clusters

    # for z in range(3): For visualizing cluster characteristics
    #    plotcluster(z)

    # average clusters within a ppm window of their suspected proteins
    sum_df = cluster_df(testq, ppm=50)

    #check mz_error kde plotting flags
    if original_mz_kde_path is not None and adjusted_mz_kde_path is not None:
        # generate plot of KDE before ppm correction
        kde_plot(sum_df, args.original_mz_kde_path)

    if calibration_outpath is not None:
        # apply polyfit mz calibration

        calib_dict = gen_mz_error_calib_output(
            testq=testq,
            calib_pk_fpath=calibration_outpath,
            polyfit_degree=polyfit_deg,
            ppm_tol=ppm_tolerance,
            int_tol=intensity_tolerance,
            cluster_corr_tol=cluster_corr_tolerance
        )
        testq["mz_mono_fix"] = apply_polyfit_cal_mz(
            polyfit_coeffs=calib_dict["polyfit_coeffs"], mz=df["mz_mono"]
        )
        testq["mz_mono_fix_round"] = np.round(testq["mz_mono_fix"].values, 3)

    else:

        # this is what is initially implemented for mz correction

        # identify major peak of abs_ppm_error clusters, apply correction to all monoisotopic mz values
        offset, offset_peak_width = find_offset(sum_df)
        if offset > 0:
            testq["mz_mono_fix"] = [
                x * (1000000 - offset) / (1000000) for x in df["mz_mono"]
            ]
            testq["mz_mono_fix_round"] = np.round(testq["mz_mono_fix"].values, 3)
        else:
            testq["mz_mono_fix"] = [
                x * (1000000 + offset) / (1000000) for x in df["mz_mono"]
            ]
            testq["mz_mono_fix_round"] = np.round(testq["mz_mono_fix"].values, 3)

        ppm_refilter = math.ceil(offset_peak_width / 2)


    # re-cluster on the adjusted MZ, same weights
    apply_cluster_weights(testq, dt_weight=5.0, rt_weight=0.6, mz_weight=0.006)

    db = DBSCAN()
    db.fit(testq[["cluster_im", "cluster_RT", "cluster_mz", "charge"]])
    clusters = db.fit_predict(testq[["cluster_im", "cluster_RT", "cluster_mz", "charge"]])
    testq["cluster"] = clusters

    # re-average clusters to single lines, check for duplicate RTs, save sum_df to outfile
    sum_df = cluster_df(testq, ppm=ppm_refilter, adjusted=True)

    if original_mz_kde_path is not None and adjusted_mz_kde_path is not None:
        # plot adjusted_mz KDE
        kde_plot(sum_df, adjusted_mz_kde_path)

    # check for duplicate RT-groups THIS MAY BE USELESS TODO
    no_duplicates, hits = find_rt_duplicates(sum_df)
    print("No rt Duplicates: " + str(no_duplicates))
    if not no_duplicates:
        print("DUPLICATES: " + hits)

    # send sum_df to main output
    sum_df.to_csv(intermediate_out_path, index=False)

if __name__ == "__main__":

    # set expected command line arguments

    # positional arguments
    parser = argparse.ArgumentParser(description="Reads an imtbx .peaks.isotopes file and creates an intermediate list of identified charged species to be used by make_library_master_list.py")
    parser.add_argument("isotopes_path", help="path/to/.peaks.isotopes file from undeuterated mzml")
    parser.add_argument("names_and_seqs_path", help="path/to/.csv with names and sequences of library proteins")
    parser.add_argument("intermediate_out_path", help="path/to/_intermediate.csv main output file")

    # optional arguments
    parser.add_argument("-p", "--plot", help="/path/to/directory/ to save original and adjusted mz-error kde plots, use instead of -o and -a")
    parser.add_argument("-o", "--original_mz_kde_path", help="/path/to/file to save original mz-error kde plots, use with -a")
    parser.add_argument("-a", "--adjusted_mz_kde_path", help="/path/to/file to save adjusted mz-error kde plots, use with -o")
    parser.add_argument("-c", "--calibration_outpath", help="/path/to/file for polyfit-calibration output, determines use of polyfit calibration") 
    parser.add_argument("-d", "--polyfit_deg", help="degree of polynomial curve to fit to mz data for non-linear correction", type=int, default=1)
    parser.add_argument("-t", "--ppm_tolerance", help="ppm error tolerance of observed to expected mz, defualt 50 ppm", type=float, default=50)
    parser.add_argument("-i", "--intensity_tolerance", help="minimum intensity to consider cluster, default 10E4", type=float, default=10000)
    parser.add_argument("-r", "--cluster_corr_tolerance", help="minimum correlation between isotope clusters to consider them redundant, default 0.99", type=float, default=0.99)
    parser.add_argument("-f", "--ppm_refilter", help="ppm error tolerance for post-mz-adjustment clusters, default 10 ppm", type=float, default=10)
    
    # parse given arguments
    args = parser.parse_args()    
    
    # check for any plotting argument
    if args.plot is not None or args.original_mz_kde_path is not None or args.adjusted_mz_kde_path is not None:
        # make explicit filenames if directory given
        if args.p is not None:
            args.o = args.p+"original_mz_kde_path.pdf"
            args.a = args.p+"adjusted_mz_kde_path.pdf"
        else:
            # require both explicit filenames
            if args.o is None or args.a is None:
                parser.print_help()
                print("Plotting with explicit paths requires both -o and -a to be set")
                sys.exit()

            # continue, we have -o and -a

    main(args.isotopes_path, args.names_and_seqs_path, args.intermediate_out_path, plot=args.plot, original_mz_kde_path=args.original_mz_kde_path, adjusted_mz_kde_path=args.adjusted_mz_kde_path, calibration_outpath=args.calibration_outpath, polyfit_deg=args.polyfit_deg, ppm_tolerance=args.ppm_tolerance, intensity_tolerance=args.intensity_tolerance, cluster_corr_tolerance=args.cluster_corr_tolerance, ppm_refilter=args.ppm_refilter)
