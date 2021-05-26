import os
import sys
import gzip
import ipdb
import pymzml
import argparse
import collections
import numpy as np
import pandas as pd


def main(mzml_path, return_flag=None, out_path=None, mzml_sum_outpath=None):
    """Generate LC Chromatogram by summing ionic current over IMS and m/Z dimensions.

    Args:
        mzml_path (string): path/to/file.mzml to be read into .tic.
        return_flag: option to return main output in python, for notebook context.
        out_path (string): option to save main output, path/to/file.tic.
        mzml_sum_outpath (string): option to save sum of mzml intensity, used in normalization between mzmls. 

    Returns:
        out_dict (dictionary): containing the below key/value pairs.
            ms1_ims_tic (np_array): LC Chromatogram as 2D numpy ndarray. Contains sum of ionic current for LC-RT and m/Z bins. 
            mzml_sum (float): sum of all mzml MS intensity.
    
    """
    drift_times = []
    scan_times = []

    # Uses mzML string pattern to find ims drift and m/z scan times
    lines = open(mzml_path, "rt").readlines()
    for line in lines:
        if ('<cvParam cvRef="MS" accession="MS:1002476" name="ion mobility drift time" value'
                in line):
            dt = line.split('value="')[1].split('"')[0]  # replace('"/>',''))
            drift_times.append(float(dt))

    for line in lines:
        if ('<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value='
                in line):
            st = line.split('value="')[1].split('"')[0]  # replace('"/>',''))
            scan_times.append(float(st))

    # drift_times = np.array(drift_times)
    # scan_times = np.array(scan_times)
    # scan_numbers=np.arange(0,len(scan_times))

    run = pymzml.run.Reader(mzml_path)
    # These are hardcoded values controlling the density of sampling - TODO: Arguments?
    mz_bins = 70
    lims = np.arange(
        600, 2020, 20
    )

    # print ((len(set(drift_times)),len(set(scan_times))))
    ms1_ims_tic = np.zeros(
        (len(set(drift_times)) * mz_bins, len(set(scan_times))), np.int)
    print(np.shape(ms1_ims_tic))

    id_appearance_count = collections.Counter()
    rtIndex = 0
    for spectrum in run:
        if spectrum["id"] == "TIC":
            continue
        spec_id = int(spectrum["id"] - 1)
        id_appearance_count[spec_id] += 1
        if id_appearance_count[spec_id] == 1:
            ims_bin = (
                spec_id % 200
            )  # Waters synapt-G2 has 200 IMS bins for each LC timepoint, TODO - make main argument with default and config variable
            specpeaks = np.array(spectrum.peaks("raw")).T
            if len(specpeaks) > 0:
                for mz_bin in range(mz_bins):
                    ms1_ims_tic[(ims_bin * mz_bins) + mz_bin, rtIndex] = int(
                        sum(specpeaks[1][(specpeaks[0] > lims[mz_bin]) &
                                         (specpeaks[0] < lims[mz_bin + 1])]))

            if ims_bin == 199:
                rtIndex += 1

                if rtIndex % 20 == 0:
                    print(rtIndex)
    mzml_sum = np.sum(ms1_ims_tic)

    if out_path is not None:
        np.savetxt(out_path, ms1_ims_tic, fmt="%i")

    if mzml_sum_outpath is not None:
        np.savetxt(mzml_sum_outpath, np.array(mzml_sum))

    if return_flag is not None:
        return {"tic": ms1_ims_tic, "mzml_sum": mzml_sum}


if __name__ == "__main__":

    # Set expected command line arguments.
    parser = argparse.ArgumentParser(
        description=
        "Sum of Total Ionic Current over IMS and m/Z dimensions, yielding an LC-Chromatogram"
    )
    parser.add_argument("mzml_path",
                        help="path/to/file for one timepoint .mzML")
    parser.add_argument("-o",
                        "--out_path",
                        help="path/to/file for main output .ims.mz.tic")
    parser.add_argument("-s",
                        "--mzml_sum_outpath",
                        help="path/to/file for <mzml>_sum.txt")
    args = parser.parse_args()

    main(args.mzml_path, out_path=args.out_path, mzml_sum_outpath=args.mzml_sum_outpath)
