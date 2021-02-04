import os
import sys
import gzip
import copy
import zlib
import ipdb
import time
import psutil
import pymzml
import numpy as np
import pandas as pd
import _pickle as cpickle
import pickle as pk
from collections import Counter


def load_pickle_file(pickle_fpath):
    """
    loads the pickle file (without any dependence on other classes or objects or functions)
    :param pickle_fpath: file path
    :return: pickle_object
    """
    with open(pickle_fpath, "rb") as file:
        pk_object = pk.load(file)
    return pk_object


def apply_polyfit_cal_mz(polyfit_coeffs, mz):
    """
    apply polyfit coeff to transform the mz values
    :param polyfit_coeffs: polyfit coefficients
    :param mz: mz values
    :return: transformed mz values
    """
    mz_corr = np.polyval(polyfit_coeffs, mz)
    return mz_corr


# TODO
print("mzML.gz: " + str(snakemake.input[1]))
library_info = pd.read_csv(snakemake.input[0])
mzml_gz = snakemake.input[1]
mzml = mzml_gz.split("/")[-1][:-3]

# find number of mzml-source timepoint for extracting RT #TODO THIS WILL HAVE TO BE HANDLED WHEN MOVING FROM MZML TO RAW - files in config[int] will not have same extension
mask = [False for i in snakemake.config["timepoints"]]
for i in range(len(snakemake.config["timepoints"])):
    if mzml in snakemake.config[snakemake.config["timepoints"][i]]:
        mask[i] = True
tp = mask.index(True)  # index of timepoint in config['timepoints']
mask = [False for i in snakemake.config[snakemake.config["timepoints"][tp]]]
for i in range(len(snakemake.config[snakemake.config["timepoints"][tp]])):
    if (
        snakemake.config[snakemake.config["timepoints"][tp]][i] == mzml
    ):  # find index of file within config[int(tp_in_seconds)]
        mask[i] = True
n_replicate = mask.index(True)

library_info["n"] = range(len(library_info))
library_info["Drift Time MS1"] = (
    library_info["im_mono"] / 200.0 * 13.781163434903
)  # 13.78116 is a hardcoded average IMS pulse time

ret_ubounds = (
    library_info["rt_group_mean_RT_%d_%d" % (tp, n_replicate)].values
    + snakemake.config["rt_radius"]
)
ret_lbounds = (
    library_info["rt_group_mean_RT_%d_%d" % (tp, n_replicate)].values
    - snakemake.config["rt_radius"]
)

dt_ubounds = library_info["Drift Time MS1"].values * (
    1 + snakemake.config["dt_radius_scale"]
)
dt_lbounds = library_info["Drift Time MS1"].values * (
    1 - snakemake.config["dt_radius_scale"]
)

output = []
drift_times = []
scan_times = []

lines = gzip.open(mzml_gz, "rt").readlines()
for line in lines:
    if (
        '<cvParam cvRef="MS" accession="MS:1002476" name="ion mobility drift time" value'
        in line
    ):
        dt = line.split('value="')[1].split('"')[0]  # replace('"/>',''))
        drift_times.append(float(dt))

for line in lines:
    if (
        '<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value='
        in line
    ):
        st = line.split('value="')[1].split('"')[0]  # replace('"/>',''))
        scan_times.append(float(st))

drift_times = np.array(drift_times)
scan_times = np.array(scan_times)
scan_numbers = np.arange(0, len(scan_times))

process = psutil.Process(os.getpid())
print(process.memory_info().rss)
msrun = pymzml.run.Reader(mzml_gz)
print(process.memory_info().rss)

starttime = time.time()

print(time.time() - starttime, mzml_gz)

hd_mass_diff = 1.006277
c13_mass_diff = 1.00335
isotope_totals = [
    len(seq) + snakemake.config["high_mass_margin"]
    for seq in library_info["sequence"].values
]

scan_to_lines = [[] for i in scan_times]
scans_per_line = []
output_scans = [[] for i in range(len(library_info))]

for i in range(len(library_info)):
    # print i
    # sys.stdout.flush()
    # TODO: This was a hack to limit pipeline runs to a subset, should be an organized feature
    if i in range(len(library_info)):  # sample_indices:
        keep_scans = scan_numbers[
            (drift_times >= dt_lbounds[i])
            & (drift_times <= dt_ubounds[i])
            & (scan_times <= ret_ubounds[i])
            & (scan_times >= ret_lbounds[i])
        ]
        scans_per_line.append(len(keep_scans))
        for scan in keep_scans:
            scan_to_lines[scan].append(i)
    else:
        scans_per_line.append(None)

    if i % 100 == 0:
        print(str(i) + " lines, time: " + str(time.time() - starttime))

relevant_scans = [i for i in scan_numbers if len(scan_to_lines[i]) > 0]
print("N Scans: " + str(len(relevant_scans)))


# implement polyfit calibration if True in config file
apply_polyfit_mz_calibration = snakemake.config["polyfit_calibration"]
calib_dict = load_pickle_file(snakemake.input[2])

for scan_number in relevant_scans:

    if scan_number % 1 == 0:
        print(
            scan_number,
            process.memory_info().rss / (1024 * 1024 * 1024),
            (len(library_info) - output_scans.count([])) / len(library_info),
        )

    if len(scan_to_lines[scan_number]) > 0:
        try:
            scan = msrun[scan_number]
            spectrum = np.array(scan.peaks("raw")).astype(np.float32)
            if len(spectrum) == 0:
                spectrum = scan.peaks("raw").astype(np.float32)
            spectrum = spectrum[spectrum[:, 1] > 10]
            # apply calibration to mz values
            if apply_polyfit_mz_calibration:
                spectrum[:, 0] = apply_polyfit_cal_mz(polyfit_coeffs=calib_dict["polyfit_coeffs"], mz=spectrum[:, 0])
        except:
            spectrum = np.array([[0, 0]])

    for i in scan_to_lines[scan_number]:
        # if len(output_scans[i]) == 0:
        print("Library Index: " + str(i) + " Len Output: " + str(len(output_scans[i])))
        obs_mz_values = library_info["obs_mz"].values[i]
        mz_low = obs_mz_values - (
            snakemake.config["low_mass_margin"] / library_info["charge"].values[i]
        )
        mz_high = obs_mz_values + (isotope_totals[i] / library_info["charge"].values[i])
        try:
            output_scans[i].append(
                spectrum[(mz_low < spectrum[:, 0]) & (spectrum[:, 0] < mz_high)]
            )
        except:
            print(i, output_scans[i], mz_low, mz_high)
            print(spectrum)
            print(spectrum[(mz_low < spectrum[:, 0]) & (spectrum[:, 0] < mz_high)])
            sys.exit(0)
        try:
            if len(output_scans[i]) == scans_per_line[i]:
                keep_drift_times = drift_times[
                    (drift_times >= dt_lbounds[i])
                    & (drift_times <= dt_ubounds[i])
                    & (scan_times <= ret_ubounds[i])
                    & (scan_times >= ret_lbounds[i])
                ]
                keep_scan_times = scan_times[
                    (drift_times >= dt_lbounds[i])
                    & (drift_times <= dt_ubounds[i])
                    & (scan_times <= ret_ubounds[i])
                    & (scan_times >= ret_lbounds[i])
                ]
                output = [
                    sorted(set(keep_scan_times)),
                    sorted(set(keep_drift_times)),
                    output_scans[i],
                ]
                # FIX FOR NEW FILE SCHEME TODO
                my_out = [
                    out
                    for out in snakemake.output
                    if out
                    == "resources/tensors/" + str(i) + "_" + mzml + ".gz.cpickle.zlib"
                ][0]
                print("My_out: " + str(my_out))
                with open(my_out, "wb") as file:
                    file.write(zlib.compress(cpickle.dumps(output)))
                print(
                    scan_number,
                    process.memory_info().rss / (1024 * 1024 * 1024),
                    "presave",
                )
                output_scans[i] = []
                print(
                    scan_number,
                    process.memory_info().rss / (1024 * 1024 * 1024),
                    "savedisk",
                )
        except:
            ipdb.set_trace()

    if len(scan_to_lines[scan_number]) > 0:
        cur_lengths = np.array(
            [len(output_scans[i]) for i in scan_to_lines[scan_number]]
        )
        target_lengths = np.array(
            [i for i in scan_to_lines[scan_number]]
        )  # np.array([scans_per_line[i] for i in scan_to_lines[scan_number]])
