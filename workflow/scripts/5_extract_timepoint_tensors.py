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
import gzip
import glob
import copy
import zlib
import time
import yaml
import psutil
import pymzml
import argparse
import numpy as np
import pandas as pd
import _pickle as cpickle
import pickle as pk
from collections import Counter
from collections import OrderedDict

def load_pickle_file(pickle_fpath):    
    """Loads a pickle file.

    Args:
        pickle_fpath (str): path/to/file.pickle

    Returns:
        pk_object (Python Object: unpickled object
    
    """
    with open(pickle_fpath, "rb") as file:
        pk_object = pk.load(file)
    return pk_object


def apply_polyfit_cal_mz(polyfit_coeffs, mz):
    """Apply mz calibration determined in make_master_list.py to an extracted tensor.

    Args:
        polyfit_coeffs (list of floats): polyfit coefficients
        mz (list of floats): mz values

    Returns:
        corrected_mz (Numpy NDarray): transformed mz values
    
    """
    corrected_mz = np.polyval(polyfit_coeffs, mz)
    return corrected_mz


def main(library_info_path,
         mzml_gz_path,
         timepoints_dict,
         outputs=None,
         return_flag=False,
         low_mass_margin=10,
         high_mass_margin=17,
         rt_radius=0.4,
         dt_radius_scale=0.06,
         polyfit_calibration_dict=None,
         indices=None):
    """Reads through .mzML file and extracts subtensors whose dimensions are defined in 
       library_info.json, optionally saves individual tensors or returns all as a dictionary.

    Args:
        library_info_path (str): path/to/library_info.json or checked_library_info.json
        mzml_gz_path (str): path/to/timepoint.mzML.gz
        timepoints_dict (dict): dictionary with 'timepoints' key containing list of hdx timepoints in integer seconds, which are keys mapping to lists of each timepoint's replicate .mzML filenames 
        outputs (list of strings): list of filename strings for writing extracted outputs. 
        return_flag (bool): option to return main output in python, for notebook context
        low_mass_margin (int): number of m/Z bins to extend the lower bound of extraction from base-peak m/Z, helps capture incompletely centered data and makes plots more readable
        high_mass_margin (int): number of m/Z bins to extend the upper bound of extraction from (base-peak + possible mass addition by number residues), helps capture incompletely centered data and makes plots more readable
        rt_radius (float): radius around signal center of mass to extract in LC - retention time
        dt_radius_scale (float): scale of radius around signal center of mass to extract in IMS - drift time
        polyfit_calibration_dict (dict): dictionary of mz-adjustment terms optionally calculated in make_library_master_list.py
        indices (list of ints): subset of library_info indices to extract

    Returns:
        out_dict (dict): dictionary containing every extracted tensor with library_info indices as keys

    """
    out_dict = {}
    library_info = pd.read_json(library_info_path)

    # Makes nested dictionary where rt-group-name is the outermost key and returns a dictionary 
    # mapping charge states of that rt-group to their library_info indices.
    names = list(OrderedDict.fromkeys(library_info["name"].values).keys()) # This is the Python-native version of an ordered set operation.
    name_charge_idx = {name: {charge: library_info.loc[(library_info["name"]==name) & (library_info["charge"]==charge)].index 
                            for charge in library_info.loc[library_info["name"]==name]["charge"].values} 
                            for name in names}

    mzml = mzml_gz_path.split("/")[-1][:-3] # Strip '.gz' from input filename to match config timepoint values.

    # Find number of mzml-source timepoint for extracting RT #TODO THIS WILL HAVE TO BE HANDLED WHEN MOVING FROM MZML TO RAW - files in config[int] will not have same extension.
    mask = [False for i in timepoints_dict["timepoints"]]
    for i in range(len(timepoints_dict["timepoints"])):
        if mzml in timepoints_dict[timepoints_dict["timepoints"][i]]:
            mask[i] = True
    tp = mask.index(True)  # index of timepoint in config['timepoints']
    mask = [False for i in timepoints_dict[timepoints_dict["timepoints"][tp]]]
    for i in range(len(timepoints_dict[timepoints_dict["timepoints"][tp]])):
        if (timepoints_dict[timepoints_dict["timepoints"][tp]][i] == mzml
           ):  # find index of file within config[int(tp_in_seconds)]
            mask[i] = True
    n_replicate = mask.index(True)

    library_info["n"] = range(len(library_info))

    # 13.78116 is a hardcoded average IMS pulse time TODO: This should be exposed to argument layer with default as well
    library_info["Drift Time MS1"] = (library_info["im_mono"] / 200.0 *
                                      13.781163434903)
    ret_ubounds = (library_info["rt_group_mean_RT_%d_%d" %
                                (tp, n_replicate)].values + rt_radius)
    ret_lbounds = (library_info["rt_group_mean_RT_%d_%d" %
                                (tp, n_replicate)].values - rt_radius)
    dt_ubounds = library_info["Drift Time MS1"].values * (1 + dt_radius_scale)
    dt_lbounds = library_info["Drift Time MS1"].values * (1 - dt_radius_scale)

    drift_times = []
    scan_times = []

    # TODO replace this hardcoded search string with an optional search string parameter with this value as default?
    lines = gzip.open(mzml_gz_path, "rt").readlines()
    for line in lines:
        if ('<cvParam cvRef="MS" accession="MS:1002476" name="ion mobility drift time" value'
                in line):
            dt = line.split('value="')[1].split('"')[0]  # replace('"/>',''))
            drift_times.append(float(dt))
    drift_times = np.array(drift_times)

    # Display memory use before beginning iteration.
    process = psutil.Process(os.getpid())
    print(process.memory_info().rss)
    msrun = pymzml.run.Reader(mzml_gz_path)
    print(process.memory_info().rss)
    starttime = time.time()
    print(time.time() - starttime, mzml_gz_path)

    scan_functions = []
    scan_times = []

    # Read scan info from msrun.
    for k in range(msrun.get_spectrum_count()):
        nextscan = msrun.next()
        scan_functions.append(nextscan.id_dict["function"])
        scan_times.append(nextscan.scan_time_in_minutes())
    scan_times = np.array(scan_times)
    scan_numbers = np.arange(0, len(scan_times))
    scan_functions = np.array(scan_functions)

    # Set upper m/Z bounds for each sequence
    isotope_totals = [
        len(seq) + high_mass_margin for seq in library_info["sequence"].values
    ]

    # Map scans to lines needing those scans, scan numbers needed by each line, and data for each line. 
    scan_to_lines = [[] for i in scan_times]
    scans_per_line = []
    output_scans = [[] for i in range(len(library_info))]

    # Use mapping lists to get scan numbers for each POI
    for i in range(len(library_info)):
        # check for subset indices
        if indices is not None:
            # only keep scans for relevant indices
            if i in indices:
                keep_scans = scan_numbers[(drift_times >= dt_lbounds[i]) &
                                          (drift_times <= dt_ubounds[i]) &
                                          (scan_times <= ret_ubounds[i]) &
                                          (scan_times >= ret_lbounds[i]) &
                                          (scan_functions == 1)]
                scans_per_line.append(len(keep_scans))
                for scan in keep_scans:
                    scan_to_lines[scan].append(i)
            else:
                scans_per_line.append(None)

        # extract for each line by default
        else:
            keep_scans = scan_numbers[(drift_times >= dt_lbounds[i]) &
                                      (drift_times <= dt_ubounds[i]) &
                                      (scan_times <= ret_ubounds[i]) &
                                      (scan_times >= ret_lbounds[i]) &
                                      (scan_functions == 1)]
            scans_per_line.append(len(keep_scans))
            for scan in keep_scans:
                scan_to_lines[scan].append(i)

        if i % 100 == 0:
            print(str(i) + " lines, time: " + str(time.time() - starttime))

    # Filter scans that don't need to be read.
    relevant_scans = [i for i in scan_numbers if len(scan_to_lines[i]) > 0]
    print("N Scans: " + str(len(relevant_scans)))

    # Perform polyfit calibration if adjustment dict passed.
    if polyfit_calibration_dict is not None:
        calib_dict = load_pickle_file(polyfit_calibration_dict)

    print(process.memory_info().rss)
    msrun = pymzml.run.Reader(mzml_gz_path)
    print(process.memory_info().rss)
    print(time.time() - starttime, mzml_gz_path)

    # Iterate over each scan, check if any line needs it, extract data within bounds for all lines needing data from the scan. 
    for scan_number in range(
            msrun.get_spectrum_count()):  
        scan = msrun.next()

        if scan_number in relevant_scans:

            # Print progress at interval.
            if scan_number % 1 == 0:
                print(
                    scan_number,
                    process.memory_info().rss / (1024 * 1024 * 1024),
                    (len(library_info) - output_scans.count([])) /
                    len(library_info),
                )
            spectrum = np.array(scan.peaks("raw")).astype(np.float32)
            if len(spectrum) == 0:
                spectrum = scan.peaks("raw").astype(np.float32)
            spectrum = spectrum[spectrum[:, 1] > 10]
            # Apply calibration to mz values if calibration dict passed.
            if polyfit_calibration_dict is not None:
                spectrum[:, 0] = apply_polyfit_cal_mz(
                    polyfit_coeffs=calib_dict["polyfit_coeffs"],
                    mz=spectrum[:, 0])

            # Iterate over each library_info index that needs to read the scan.
            for i in scan_to_lines[scan_number]:  
                print("Library Index: " + str(i) + " Len Output: " +
                      str(len(output_scans[i])))
                obs_mz_values = library_info["obs_mz"].values[i]
                mz_low = obs_mz_values - (low_mass_margin /
                                          library_info["charge"].values[i])
                mz_high = obs_mz_values + (isotope_totals[i] /
                                           library_info["charge"].values[i])
                try:
                    output_scans[i].append(spectrum[(mz_low < spectrum[:, 0]) &
                                                    (spectrum[:, 0] < mz_high)])
                except:
                    print("spectrum read error, scan: " + str(scan_number) +
                          " , line: " + str(i))
                    print(i, output_scans[i], mz_low, mz_high)
                    print(spectrum)
                    print(spectrum[(mz_low < spectrum[:, 0]) &
                                   (spectrum[:, 0] < mz_high)])
                    sys.exit(0)
                #try:
                # Is this the last scan the line needed? If so, save to disk.
                if len(output_scans[i]) == scans_per_line[i]:
                    my_name = library_info.iloc[i]["name"]
                    my_charge = library_info.iloc[i]["charge"]
                    my_out_test = "/" + my_name + "/" + my_name + "_" + "charge" + str(my_charge) + "_" + mzml + ".gz.cpickle.zlib"
                    print("name = " + str(my_name), flush=True)
                    print("charge = " + str(my_charge), flush=True)
                    print("my_out search string = " + my_out_test, flush=True)
                    keep_drift_times = drift_times[
                        (drift_times >= dt_lbounds[i]) &
                        (drift_times <= dt_ubounds[i]) &
                        (scan_times <= ret_ubounds[i]) &
                        (scan_times >= ret_lbounds[i]) &
                        (scan_functions == 1)]
                    keep_scan_times = scan_times[
                        (drift_times >= dt_lbounds[i]) &
                        (drift_times <= dt_ubounds[i]) &
                        (scan_times <= ret_ubounds[i]) &
                        (scan_times >= ret_lbounds[i]) &
                        (scan_functions == 1)]
                    output = [
                        sorted(set(keep_scan_times)),
                        sorted(set(keep_drift_times)),
                        output_scans[i],
                    ]

                    if return_flag:
                        out_dict[i] = output

                    # Save to file if outputs provided, match name and charge to index. 
                    if outputs is not None:
                        my_out = [
                            out for out in outputs if "/" + my_name + "/" + my_name 
                            + "_" + "charge" + str(my_charge) + "_" + mzml + ".gz.cpickle.zlib" in out
                        ][0]
                        print("my_out: " + str(my_out), flush=True)
                        with open(my_out, "wb+") as file:
                            file.write(zlib.compress(cpickle.dumps(output)))
                        print(
                            scan_number,
                            process.memory_info().rss /
                            (1024 * 1024 * 1024),
                            "presave",
                        )
                        output_scans[i] = []
                        print(
                            scan_number,
                            process.memory_info().rss /
                            (1024 * 1024 * 1024),
                            "savedisk",
                        )
                    else:
                        output_scans[i] = []  # Avoid duplication of tensors in return-only state.
                """
                except:
                    print("error in output block on scan: " + str(scan_number) +
                          " , for line: " + str(i))
                    sys.stdout.flush()
                    sys.exit(0)
                """
    if return_flag:
        return out_dict


if __name__ == "__main__":

    # If the snakemake global object is present, save expected arguments from snakemake to be passed to main().
    if "snakemake" in globals():
        polyfit_calibration_dict = None
        indices = None
        open_timepoints = yaml.load(open(snakemake.input[2], "rb").read(), Loader=yaml.Loader)
        # Check for optional arguments.
        if len(snakemake.input) > 3:
            if ".pk" in snakemake.input[3]:
                polyfit_calibration_dict = snakemake.input[3]
                if len(snakemake.input) > 4:
                    indices = pd.read_csv(snakemake.input[4])['index'].values
            else:
                indices = pd.read_csv(snakemake.input[3])['index'].values

        # Obtain dt and rt radius from config file.
        config_rt_radius = open_timepoints["rt_radius"]
        config_dt_radius_scale = open_timepoints["dt_radius_scale"]

        main(library_info_path=snakemake.input[0],
             mzml_gz_path=snakemake.input[1],
             timepoints_dict=open_timepoints,
             outputs=snakemake.output,
             rt_radius = config_rt_radius,
             dt_radius_scale = config_dt_radius_scale,
             polyfit_calibration_dict=polyfit_calibration_dict,
             indices=indices)
    else:
        # CLI context, set expected arguments with argparse module.
        parser = argparse.ArgumentParser()
        parser.add_argument("library_info_path", help="path/to/library_info.json")
        parser.add_argument("mzml_gz_path", help="path/to/file.mzML.gz")
        parser.add_argument(
            "timepoints_yaml",
            help=
            "path/to/file.yaml containing list of hdx timepoints in integer seconds which are also keys mapping to lists of each timepoint's .mzML file, can pass config/config.yaml - for Snakemake context"
        )
        parser.add_argument(
            "-u",
            "--high_mass_margin",
            default=17,
            help=
            "radius around expected rt to extend extraction window in rt-dimension")
        parser.add_argument(
            "-l",
            "--low_mass_margin",
            default=10,
            help=
            "integrated-mz-bin magnitude of margin behind the POI monoisotopic mass, to avoid signal truncation"
        )
        parser.add_argument(
            "-r",
            "--rt_radius",
            default=0.4,
            help=
            "integrated-m/z-bin magnitude of margin beyond estimated full-deuteration, to avoid signal truncation"
        )
        parser.add_argument(
            "-d",
            "--dt_radius_scale",
            default=0.06,
            help=
            "scale factor for radius around expected dt to extend extraction window in dt-dimension"
        )
        parser.add_argument(
            "-c",
            "--polyfit_calibration_dict",
            help=
            "path/to/file_mz_calib_dict.pk, provide if using polyfit mz recalibration"
        )
        parser.add_argument("-o",
                            "--outputs",
                            nargs="*",
                            help="explicit list of string outputs to be created")
        parser.add_argument(
            "-i",
            "--indices_csv",
            help="filter_passing_indices.csv with 'index' argument, subset of library_info to extract tensors for, use with -o or -t")
        parser.add_argument(
            "-t",
            "--output_directory",
            help=
            "path/to/output_dir/ to generate outputs automatically, using without -i will extract all charged species from library_info, overridden by -o"
        )
        args = parser.parse_args()
        
        # Handle implicit arguments.
        if args.outputs is None:
            if args.output_directory is None:
                parser.print_help()
                sys.exit()
            else:
                library_info = pd.read_json(args.library_info_path)
                mzml = args.mzml_gz_path.split("/")[-1][:-3]
                if args.indices_csv is not None:
                    indices = pd.read_csv(args.indices_csv)["index"].values
                    # Only make subdirs for indices to be used
                    for i in indices:
                        if not os.path.isdir(args.output_directory+str(i)+"/"):
                            os.mkdir(args.output_directory+str(i)+"/")
                    # Make subset outputs.
                    args.outputs = [
                        args.output_directory + str(i) + "/" + str(i) + "_" + mzml + ".gz.cpickle.zlib"
                        for i in args.indices
                    ]
                else:
                    # Make all subdirs if needed.
                    for i in range(len(library_info)):
                        if not os.path.isdir(args.output_directory+str(i)+"/"):
                            os.mkdir(args.output_directory+str(i)+"/")
                    # Make an output for each line in library_info. TODO: Add an indices check.
                    args.outputs = [
                        args.output_directory + str(i) + "/" + str(i) + "_" + mzml + ".gz.cpickle.zlib"
                        for i in range(len(library_info))
                    ]

        open_timepoints = yaml.load(open(args.timepoints_yaml, "rb").read(), Loader=yaml.Loader)

        # Obtain dt and rt radius from config file if not passed to argparse.
        if args.rt_radius is None:
            args.rt_radius = open_timepoints["rt_radius"]
        if args.dt_radius_scale is None:
            args.dt_radius_scale = open_timepoints["dt_radius_scale"]

        main(library_info_path=args.library_info_path,
             mzml_gz_path=args.mzml_gz_path,
             timepoints_dict=open_timepoints,
             outputs=args.outputs,
             low_mass_margin=args.low_mass_margin,
             high_mass_margin=args.high_mass_margin,
             rt_radius=args.rt_radius,
             dt_radius_scale=args.dt_radius_scale,
             polyfit_calibration_dict=args.polyfit_calibration_dict,
             indices=args.indices)
