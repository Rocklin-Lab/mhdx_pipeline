import os
import sys
import copy
import math
import argparse
import numpy as np
import pandas as pd

sys.path.append(os.getcwd() + "/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx

def optimize_paths_inputs(library_info_path, input_directory_path, rt_group_name, timepoints):
    """Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value

   """

    # Pass inputs as fxn of rt-group name. Creates _tensor() input filenames in fixed pattern, input tensor names include library_info.index and rt-group avg elution time.
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
                    )  # TODO: This may break when using .raw as first input, investigate
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

def main(library_info_path=None, all_tensor_input_paths=None, input_directory_path=None, rt_group_name=None, timepoints=None, old_data_dir=None, html_plot_out_path=None, winner_out_path=None, runner_out_path=None, undeut_ground_out_path=None, winner_scores_out_path=None, rtdt_com_cvs_out_path=None):
    """Summary or Description of the Function

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value

   """

    # open library_info
    library_info = pd.read_csv(library_info_path)

    # order files, pooling all replicates and charges by timepoint
    name = snakemake.wildcards["name"]
    atc = []
    for tp in timepoints["timepoints"]:
        tp_buf = []
        for fn in timepoints[tp]:
            for file in all_tensor_input_paths:
                if fn in file:
                    ics = hx.limit_read(file)  # expects list of ics
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
        hx.limit_write(p1.winner, winner_out_path)
    if runner_out_path is not None:
        hx.limit_write(p1.runners, runner_out_path)
    if undeut_ground_out_path is not None:
        hx.limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products], undeut_ground_out_path)
    if winner_scores_out_path is not None:
        hx.limit_write(p1.winner_scores, winner_scores_out_path)
    if rtdt_com_cvs_out_path is not None:
        hx.limit_write([p1.rt_com_cv, p1.dt_com_cv], rtdt_com_cvs_out_path)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    # inputs
    parser.add_argument("library_info_path", help="path/to/library_info.csv")
    parser.add_argument("-i", "--all_tensor_input_paths", nargs='*', help="structured 2D list of extracted tensor file paths by HDX timepoint, include all charges and all timepoint replicates in each timepoint list.")
    parser.add_argument("-d", "--input_directory_path", help="path/to/directory to search for relevant files if assembling filenames automatically, requires -n")
    parser.add_argument("-n", "--rt_group_name", help="rt-group name to use for generating relevant tensor files, requires -d")
    parser.add_argument("-t", "--timepoints", help="dictionary with 'timepoints' containing hdx times in seconds, and a key for each timepoint corresponding to a list of timepoint mzml filenames. Can pass opened snakemake.config object")
    parser.add_argument("-g", "--old_data_dir", help="directory containing Gabe's pickled output files, using this option prints old data on plots")
    # outputs
    parser.add_argument("-o", "--html_plot_out_path", help="path/to/file for .html plot of results")
    parser.add_argument("-w", "--winner_out_path", help="path/to/file to save winning IsotopeCluster objects")
    parser.add_argument("-r", "--runner_out_path", help="path/to/file to save runner-up IsotopeClusters")
    parser.add_argument("-p", "--undeut_ground_out_path", help="path/to/file to save selected highest-confidence undeuterated IsotopeClusters")
    parser.add_argument("-s", "--winner_scores_out_path", help="path/to/file to save winning path IC scores")
    parser.add_argument("-c", "--rtdt_com_cvs_out_path", help="path/to/file to save rt/dt error measurement")

    args = parser.parse_args()

    #check for explicit inputs
    if args.all_tensor_input_paths is None:
        if args.input_directory_path is not None and args.rt_group_name is not None:
            args.all_tensor_input_paths = optimize_paths_inputs(args.library_info_path, args.input_directory_path, args.rt_group_name, args.timepoints)
        else:
            parser.print_help()
            sys.exit()

    main(library_info_path=args.library_info_path, all_tensor_input_paths=args.all_tensor_input_paths, input_directory_path=args.input_directory_path, rt_group_name=args.rt_group_name, timepoints=args.timepoints, old_data_dir=args.old_data_dir, html_plot_out_path=args.html_plot_out_path, winner_out_path=args.winner_out_path, runner_out_path=args.runner_out_path, undeut_ground_out_path=args.undeut_ground_out_path, winner_scores_out_path=args.winner_scores_out_path, rtdt_com_cvs_out_path=args.rtdt_com_cvs_out_path)
    