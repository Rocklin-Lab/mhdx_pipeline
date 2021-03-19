import os
import sys
import gzip
import glob
import copy
import zlib
import ipdb
import time
import psutil
import pymzml
import argparse
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


def main(library_info_path, mzml_gz_path, timepoints, outputs, low_mass_margin=10, high_mass_margin=17, rt_radius=0.4, dt_radius_scale=0.06, polyfit_calibration_dict=None, indices=None):
    library_info = pd.read_csv(library_info_path)
    mzml = mzml_gz_path.split("/")[-1][:-3]

    # find number of mzml-source timepoint for extracting RT #TODO THIS WILL HAVE TO BE HANDLED WHEN MOVING FROM MZML TO RAW - files in config[int] will not have same extension
    mask = [False for i in timepoints["timepoints"]]
    for i in range(len(timepoints["timepoints"])):
        if mzml in timepoints[timepoints["timepoints"][i]]:
            mask[i] = True
    tp = mask.index(True)  # index of timepoint in config['timepoints']
    mask = [False for i in timepoints[timepoints["timepoints"][tp]]]
    for i in range(len(timepoints[timepoints["timepoints"][tp]])):
        if (
            timepoints[timepoints["timepoints"][tp]][i] == mzml
        ):  # find index of file within config[int(tp_in_seconds)]
            mask[i] = True
    n_replicate = mask.index(True)

    library_info["n"] = range(len(library_info))
    library_info["Drift Time MS1"] = (
        library_info["im_mono"] / 200.0 * 13.781163434903
    )  # 13.78116 is a hardcoded average IMS pulse time TODO: This should be exposed to argument layer with default as well

    ret_ubounds = (
        library_info["rt_group_mean_RT_%d_%d" % (tp, n_replicate)].values
        + rt_radius
    )
    ret_lbounds = (
        library_info["rt_group_mean_RT_%d_%d" % (tp, n_replicate)].values
        - rt_radius
    )

    dt_ubounds = library_info["Drift Time MS1"].values * (
        1 + dt_radius_scale
    )
    dt_lbounds = library_info["Drift Time MS1"].values * (
        1 - dt_radius_scale
    )

    output = []
    drift_times = []
    scan_times = []

    lines = gzip.open(mzml_gz_path, "rt").readlines()
    for line in lines:
        if (
            '<cvParam cvRef="MS" accession="MS:1002476" name="ion mobility drift time" value'
            in line
        ):
            dt = line.split('value="')[1].split('"')[0]  # replace('"/>',''))
            drift_times.append(float(dt))
    drift_times = np.array(drift_times)

    #for line in lines:
    #    if (
    #        '<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value='
    #        in line
    #    ):
    #        st = line.split('value="')[1].split('"')[0]  # replace('"/>',''))
    #        scan_times.append(float(st))
    #scan_times = np.array(scan_times)
    #scan_numbers = np.arange(0, len(scan_times))

    process = psutil.Process(os.getpid())
    print(process.memory_info().rss)
    msrun = pymzml.run.Reader(mzml_gz_path)
    print(process.memory_info().rss)
    starttime = time.time()
    print(time.time() - starttime, mzml_gz_path)

    scan_functions = []
    scan_times = []


    for k in range(msrun.get_spectrum_count()):
        nextscan = msrun.next()
        scan_functions.append(nextscan.id_dict['function'])
        scan_times.append(nextscan.scan_time_in_minutes())

    scan_times = np.array(scan_times)
    scan_numbers = np.arange(0, len(scan_times))
    scan_functions = np.array(scan_functions)



    hd_mass_diff = 1.006277
    c13_mass_diff = 1.00335
    isotope_totals = [
        len(seq) + high_mass_margin
        for seq in library_info["sequence"].values
    ]

    scan_to_lines = [[] for i in scan_times]
    scans_per_line = []
    output_scans = [[] for i in range(len(library_info))]

    for i in range(len(library_info)):
        # Check for subset indices
        if indices is not None:
            # only keep scans for relevant indices 
            if i in indices:
                keep_scans = scan_numbers[
                    (drift_times >= dt_lbounds[i])
                    & (drift_times <= dt_ubounds[i])
                    & (scan_times <= ret_ubounds[i])
                    & (scan_times >= ret_lbounds[i])
                    & (scan_functions == 1)
                ]
                scans_per_line.append(len(keep_scans))
                for scan in keep_scans:
                    scan_to_lines[scan].append(i)
            else:
                scans_per_line.append(None)
        
        # extract for each line by default
        else:
            keep_scans = scan_numbers[
                (drift_times >= dt_lbounds[i])
                & (drift_times <= dt_ubounds[i])
                & (scan_times <= ret_ubounds[i])
                & (scan_times >= ret_lbounds[i])
                & (scan_functions == 1)
            ]
            scans_per_line.append(len(keep_scans))
            for scan in keep_scans:
                scan_to_lines[scan].append(i)

        if i % 100 == 0:
            print(str(i) + " lines, time: " + str(time.time() - starttime))

    relevant_scans = [i for i in scan_numbers if len(scan_to_lines[i]) > 0]
    print("N Scans: " + str(len(relevant_scans)))


    # implement polyfit calibration if True in config file
    if polyfit_calibration_dict is not None: 
        calib_dict = load_pickle_file(polyfit_calibration_dict)

    print(process.memory_info().rss)
    msrun = pymzml.run.Reader(mzml_gz_path)
    print(process.memory_info().rss)
    print(time.time() - starttime, mzml_gz_path)

    for scan_number in range(msrun.get_spectrum_count()):
        scan = msrun.next()
        if scan_number in relevant_scans:
        
            if scan_number % 1 == 0:
                print(
                    scan_number,
                    process.memory_info().rss / (1024 * 1024 * 1024),
                    (len(library_info) - output_scans.count([])) / len(library_info),
                )

            if len(scan_to_lines[scan_number]) > 0:
                try:
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
                    low_mass_margin / library_info["charge"].values[i]
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
                            & (scan_functions == 1)
                        ]
                        keep_scan_times = scan_times[
                            (drift_times >= dt_lbounds[i])
                            & (drift_times <= dt_ubounds[i])
                            & (scan_times <= ret_ubounds[i])
                            & (scan_times >= ret_lbounds[i])
                            & (scan_functions == 1)
                        ]
                        output = [
                            sorted(set(keep_scan_times)),
                            sorted(set(keep_drift_times)),
                            output_scans[i],
                        ]
                        # FIX FOR NEW FILE SCHEME TODO
                        my_out = [
                            out
                            for out in outputs
                            if "/" + str(i) + "_" + mzml + ".gz.cpickle.zlib" in out
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

if __name__ == "__main__":

    # set expected command line arguments
    parser = argparse.ArgumentParser()
    # inputs
    parser.add_argument("library_info_path", help="path/to/library_info.csv")
    parser.add_argument("mzml_gz_path", help="path/to/file.mzML.gz")
    parser.add_argument("timepoints", help="dictionary with 'timepoints' containing hdx times in seconds, and a key for each timepoint corresponding to a list of timepoint mzml filenames. Can pass opened snakemake.config object")
    parser.add_argument("-u", "--high_mass_margin", default=17, help="radius around expected rt to extend extraction window in rt-dimension")
    parser.add_argument("-l", "--low_mass_margin", default=10, help="integrated-mz-bin magnitude of margin behind the POI monoisotopic mass, to avoid signal truncation")
    parser.add_argument("-r", "--rt_radius", defualt=0.4, help="integrated-m/z-bin magnitude of margin beyond estimated full-deuteration, to avoid signal truncation")
    parser.add_argument("-d", "--dt_radius_scale", default=0.06, help="scale factor for radius around expected dt to extend extraction window in dt-dimension")
    parser.add_argument("-c", "--polyfit_calibration_dict", help="path/to/file_mz_calib_dict.pk, provide if using polyfit mz recalibration")
    # outputs
    parser.add_argument("-o", "--outputs", nargs='*', help="explicit list of string outputs to be created")
    parser.add_argument("-i", "--indices", nargs='*', type=int, help="subset of library_info to extract tensors for, use with -o or -t")
    parser.add_argument("-t", "--output_directory", help="path/to/output_dir/ to generate outputs automatically, using without -i will extract all charged species from library_info, overridden by -o")
    # parse given arguments
    args = parser.parse_args()

    # generate explicit output paths if not given
    if args.outputs is None:
        if args.output_directory is None:
            parser.print_help()
            sys.exit()
        else:
            if args.indices is not None:
                library_info = pd.read_csv(args.library_info_path)
                mzml = args.mzml_gz_path.split("/")[-1][:-3]
                args.outputs = [args.output_directory+str(i)+"_"+mzml for i in args.indices]
                
    main(library_info_path=args.library_info_path, mzml_gz_path=args.mzml_gz_path, timepoints=args.timepoints, outputs=args.outputs, low_mass_margin=args.low_mass_margin, high_mass_margin=args.high_mass_margin, rt_radius=args.rt_radius, dt_radius_scale=args.dt_radius_scale, polyfit_calibration_dict=args.polyfit_calibration_dict, indices=args.indices)






