import os
import sys
import copy
import math
import argparse
import numpy as np
import pandas as pd

sys.path.append(os.getcwd() + "/workflow/scripts/auxiliary/")
from HDX_LIMIT import limit_read, limit_write, PathOptimizer

def optimize_paths_inputs(library_info_path, input_directory_path, rt_group_name, timepoints):
    """Generate explicit PathOptimizer input paths for one rt_group

    Parameters:
    library_info_path (str): path/to/library_info.csv
    input_directory_path (str): /path/to/dir/ to prepend to each outpath
    rt_group_name (str): value from 'name' column of library_info
    timepoints (dict): dictionary containing list of hdx timepoints in seconds, where each timepoint is also an integer key corresponding to that timepoint's .mzML filenames

    Returns:
    name_inputs (list of strings): flat list of all IsotopeCluster inputs to PathOptimizer

   """

    name_inputs = []
    library_info = pd.read_csv(library_info_path)
    idxs = library_info.index[library_info["name"] == name].tolist()
    for key in timepoints["timepoints"]:
        if len(timepoints[key]) > 1:
            for file in timepoints[key]:
                for idx in idxs:
                    name_inputs.append(
                        input_directory_path
                        + str(idx)
                        + "_"
                        + file
                        + ".gz.cpickle.zlib"
                    )  # TODO: This won't work if we go to .mzML .RAW interoperability, 
        else:
            file = timepoints[key][0]
            for idx in idxs:
                name_inputs.append(
                    input_directory_path
                    + str(idx)
                    + "_"
                    + file
                    + ".gz.cpickle.zlib"
                )

    return name_inputs

def main(library_info_path, all_tensor_input_paths, timepoints, return_flag=False, rt_group_name=None, old_data_dir=None, html_plot_out_path=None, winner_out_path=None, runner_out_path=None, undeut_ground_out_path=None, winner_scores_out_path=None, rtdt_com_cvs_out_path=None):
    """Uses PathOptimzier class to generate best-estimate of hdx-timeseries of IsotopeClusters for a given library protein

    Parameters:
    library_info_path (str): path/to/library_info.csv
    all_tensor_input_paths (list of strings): list of paths/to/file for all lists of IsotopeClusters from generate_tensor_ics.py
    timepoints (dict): dictionary with 'timepoints' key containing list of hdx timepoints in integer seconds, which are keys mapping to lists of each timepoint's replicate .mzML filenames 
    return_flag: option to return main output in python, for notebook context
    rt_group_name (str): library_info['name'] value
    old_data_dir (str): path/to/dir to provide comparison to GJR formatted results
    html_plot_out_path (str): path/to/file.html for interactive bokeh plot
    winner_out_path (str): path/to/file for winning path
    runner_out_path (str): path/to/file for top n_runners paths
    undeut_ground_out_path (str): path/to/file for undeuterated ground-truth IsotopeClusters
    winner_scores_out_path (str): path/to/file for winner path score values
    rtdt_com_cvs_out_path (str): path/to/file for rt and dt correlation values

    Returns:
    out_dict (dict): dictionary containing 'path_optimizer': PathOptimizer instance

   """

   out_dict = {}

    # open library_info
    library_info = pd.read_csv(library_info_path)

    # order files, pooling all replicates and charges by timepoint
    name = library_info_path.split('/')[-1].split('.')[0]
    atc = []
    for tp in timepoints["timepoints"]:
        tp_buf = []
        for fn in timepoints[tp]:
            for file in all_tensor_input_paths:
                if fn in file:
                    ics = limit_read(file)  # expects list of ics
                    for ic in ics:
                        tp_buf.append(ic)

        atc.append(tp_buf)


    p1 = hx.PathOptimizer(
        name,
        atc,
        library_info,
        timepoints=timepoints["timepoints"],
        n_undeut_runs=len(timepoints[0]),
        old_data_dir=old_data_dir,
    )
    
    p1.optimize_paths()
    
    # Save by options
    if html_plot_out_path is not None:
        p1.bokeh_plot(html_plot_out_path)
    if winner_out_path is not None:
        limit_write(p1.winner, winner_out_path)
    if runner_out_path is not None:
        limit_write(p1.runners, runner_out_path)
    if undeut_ground_out_path is not None:
        limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products], undeut_ground_out_path)
    if winner_scores_out_path is not None:
        limit_write(p1.winner_scores, winner_scores_out_path)
    if rtdt_com_cvs_out_path is not None:
        limit_write([p1.rt_com_cv, p1.dt_com_cv], rtdt_com_cvs_out_path)

    if return_flag:
        out_dict['path_optimizer'] = p1
        return out_dict

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Generate a best-estimate HDX-timeseries of IsotopeClusters for a given library protein")
    # inputs
    parser.add_argument("library_info_path", help="path/to/library_info.csv")
    parser.add_argument("timepoints_yaml", help="path/to/file.yaml containing list of hdx timepoints in integer seconds which are also keys mapping to lists of each timepoint's .mzML file, can pass config/config.yaml - for Snakemake context")
    parser.add_argument("-i", "--all_tensor_input_paths", nargs='*', help="structured 2D list of extracted tensor file paths by HDX timepoint, include all charges and all timepoint replicates in each timepoint list.")
    parser.add_argument("-d", "--input_directory_path", help="path/to/directory to search for relevant files if assembling filenames automatically, requires --rt_group_name")
    parser.add_argument("-n", "--rt_group_name", help="rt-group name to use for generating relevant tensor files, requires --input_directory_path")
    parser.add_argument("-g", "--old_data_dir", help="directory containing Gabe's pickled output files, using this option prints old data on plots")
    # outputs
    parser.add_argument("-o", "--html_plot_out_path", help="path/to/file for .html plot of results")
    parser.add_argument("-w", "--winner_out_path", help="path/to/file to save winning IsotopeCluster objects")
    parser.add_argument("-r", "--runner_out_path", help="path/to/file to save runner-up IsotopeClusters")
    parser.add_argument("-p", "--undeut_ground_out_path", help="path/to/file to save selected highest-confidence undeuterated IsotopeClusters")
    parser.add_argument("-s", "--winner_scores_out_path", help="path/to/file to save winning path IC scores")
    parser.add_argument("-c", "--rtdt_com_cvs_out_path", help="path/to/file to save rt/dt error measurement")

    args = parser.parse_args()

    #open .yaml into dict for main()
    open_timepoints = yaml.load(open(args.timepoints_yaml, 'rb').read())

    #check for explicit inputs
    if args.all_tensor_input_paths is None:
        if args.input_directory_path is not None and args.rt_group_name is not None:
            args.all_tensor_input_paths = optimize_paths_inputs(args.library_info_path, args.input_directory_path, args.rt_group_name, args.timepoints)
        else:
            parser.print_help()
            sys.exit()

    main(library_info_path=args.library_info_path, all_tensor_input_paths=args.all_tensor_input_paths, timepoints=open_timepoints, rt_group_name=args.rt_group_name, old_data_dir=args.old_data_dir, html_plot_out_path=args.html_plot_out_path, winner_out_path=args.winner_out_path, runner_out_path=args.runner_out_path, undeut_ground_out_path=args.undeut_ground_out_path, winner_scores_out_path=args.winner_scores_out_path, rtdt_com_cvs_out_path=args.rtdt_com_cvs_out_path)
    