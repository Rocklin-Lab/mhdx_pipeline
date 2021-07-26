import os, shutil
import sys
import copy
import math
import argparse
import numpy as np
import pandas as pd
import yaml

sys.path.append(os.getcwd() + "/workflow/scripts/")
from HDX_LIMIT.io import limit_read, limit_write
from HDX_LIMIT.processing import PathOptimizer
from HDX_LIMIT.gjr_plot import plot_gjr_


def optimize_paths_inputs(library_info_path, input_directory_path,
                          rt_group_name, timepoints):
    """Generate explicit PathOptimizer input paths for one rt_group.
    Args:
        library_info_path (str): path/to/checked_library_info.json
        input_directory_path (str): /path/to/dir/ to prepend to each outpath
        rt_group_name (str): value from 'name' column of library_info
        timepoints (dict): dictionary containing list of hdx timepoints in seconds, where each timepoint is also an integer key corresponding to that timepoint's .mzML filenames
    Returns:
        name_inputs (list of strings): flat list of all IsotopeCluster inputs to PathOptimizer
    """
    name_inputs = []
    library_info = pd.read_json(library_info_path)
    idxs = library_info.index[library_info["name"] == rt_group_name].tolist()
    for key in timepoints["timepoints"]:
        if len(timepoints[key]) > 1:
            for file in timepoints[key]:
                for idx in idxs:
                    name_inputs.append(
                        input_directory_path + str(idx) + "_" + file +
                        ".gz.cpickle.zlib"
                    )  # TODO: This won't work if we go to .mzML .RAW interoperability,
        else:
            file = timepoints[key][0]
            for idx in idxs:
                name_inputs.append(input_directory_path + str(idx) + "_" +
                                   file + ".gz.cpickle.zlib")

    return name_inputs


def gen_correlate_matrix(list_of_arrays):
    """
    generate a correlation matrix
    :param list_of_arrays: list of arrays
    :return: correlation matrix
    """

    corr_matrix = np.zeros((len(list_of_arrays), len(list_of_arrays)))
    for ind1, arr1 in enumerate(list_of_arrays):
        for ind2, arr2 in enumerate(list_of_arrays):
            corr_matrix[ind1, ind2] = max(np.correlate(arr1/np.linalg.norm(arr1), arr2/np.linalg.norm(arr2)))

    return corr_matrix


def main(library_info_path,
         all_ic_input_paths,
         timepoints,
         return_flag=False,
         rt_group_name=None,
         old_data_dir=None,
         path_plot_out_path=None,
         html_plot_out_path=None,
         winner_out_path=None,
         runner_out_path=None,
         undeut_ground_out_path=None,
         winner_scores_out_path=None,
         rtdt_com_cvs_out_path=None,
         out_monobody_path=False):
    """Uses PathOptimzier class to generate best-estimate hdx-timeseries of IsotopeClusters for a given library protein.
    Args:
        library_info_path (str): path/to/checked_library_info.json
        all_ic_input_paths (list of strings): list of paths/to/files.cpickle.zlib for all lists of IsotopeClusters from generate_tensor_ics.py
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
        out_dict (dict): dictionary containing 'path_optimizer' key, and corresponding PathOptimizer object 
    """
    out_dict = {}

    # open library_info
    library_info = pd.read_json(library_info_path)

    if rt_group_name is None:
        name = library_info.iloc[int(
            all_ic_input_paths[0].split("/")[-1].split("_")[0])]["name"]
    else:
        name = rt_group_name

    # order files, pooling all replicates and charges by timepoint
    atc = []
    for tp in timepoints["timepoints"]:
        tp_buf = []
        for fn in timepoints[tp]:
            for file in all_ic_input_paths:
                if fn.split(
                        ".")[0] in file:  # only match filename without .mzML
                    ics = limit_read(file)  # expects list of ics
                    for ic in ics:
                        tp_buf.append(ic)

        atc.append(tp_buf)

    # Generate nearest_neighbor_correlation for each ic
    # Maybe atc should be saved as a gz.cpickle.zlib file?
    for ics in atc:
        all_baseline_integrated_mz = []
        all_rts = []
        charge_list = []
        for ic in ics:
            all_baseline_integrated_mz.append(ic.baseline_integrated_mz)
            all_rts.append(ic.rts)
            charge_list.append(ic.charge_states[0])

        charge_list = np.array(charge_list)
        mz_corrmat = gen_correlate_matrix(all_baseline_integrated_mz)
        rt_corrmat = gen_correlate_matrix(all_rts)
        minimum_corrmat = np.minimum(mz_corrmat, rt_corrmat)
        for column, ic in enumerate(ics):
            min_corr_list = minimum_corrmat[column][charge_list != ic.charge_states[0]]
            if len(min_corr_list) != 0:
                ic.nearest_neighbor_correlation = max(min_corr_list)
            else:
                ic.nearest_neighbor_correlation = 0

    p1 = PathOptimizer(
        name,
        atc,
        library_info,
        timepoints=timepoints["timepoints"],
        n_undeut_runs=len(timepoints[0]),
        old_data_dir=old_data_dir,
    )

    if out_monobody_path:

        # Create folders for monobody score function
        if not os.path.isdir('/'.join(winner_out_path.split('/')[:-1])):
            os.mkdir('/'.join(winner_out_path.split('/')[:-1]))
        if not os.path.isdir('resources/ic_time_series/monobody'):
            os.mkdir('resources/ic_time_series/monobody')
        if not os.path.isdir('results/plots/ic_time_series/winner_plots/monobody'):
            os.mkdir('results/plots/ic_time_series/winner_plots/monobody')

        # Generate best paths for monobody score function
        p1.optimize_paths_mono()

        if winner_out_path is not None:
            limit_write(p1.winner, winner_out_path)
        if runner_out_path is not None:
            limit_write(p1.runners, runner_out_path)
        if undeut_ground_out_path is not None:
            limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products],
                        undeut_ground_out_path)
        if winner_scores_out_path is not None:
            limit_write(p1.winner_scores, winner_scores_out_path)
        if rtdt_com_cvs_out_path is not None:
            limit_write([p1.rt_com_cv, p1.dt_com_cv], rtdt_com_cvs_out_path)
        if path_plot_out_path is not None:
            undeut_grounds = [p1.undeut_grounds, p1.undeut_ground_dot_products]
            plot_gjr_(winner=p1.winner,
                      undeut_grounds=undeut_grounds,
                      output_path=path_plot_out_path,
                      prefix=name)

        # Transfer folder and files to monobody folder
        if not os.path.isfile('resources/ic_time_series/monobody/' + '/'.join(winner_out_path.split('/')[-2:])):
            shutil.move('/'.join(winner_out_path.split('/')[:-1]), 'resources/ic_time_series/monobody/')
        if not os.path.isfile('resources/ic_time_series/monobody/' + '/'.join(winner_out_path.split('/')[-2:])):
            shutil.move(path_plot_out_path, 'results/plots/ic_time_series/winner_plots/monobody/')

    # Create folders and move files for multibody score terms
    # Somewhat redundant to snakemake that already create those folders
    if not os.path.isdir('/'.join(winner_out_path.split('/')[:-1])):
        os.mkdir('/'.join(winner_out_path.split('/')[:-1]))
    if not os.path.isdir('resources/ic_time_series/multibody'):
        os.mkdir('resources/ic_time_series/multibody/')
    if not os.path.isdir('results/plots/ic_time_series/winner_plots/multibody'):
        os.mkdir('results/plots/ic_time_series/winner_plots/multibody')

    # Generate best paths for multibody score function
    p1.optimize_paths_multi()

    # write outputs
    # if html_plot_out_path is not None:
    #     p1.bokeh_plot(html_plot_out_path)
    if winner_out_path is not None:
        limit_write(p1.winner, winner_out_path)
    if runner_out_path is not None:
        limit_write(p1.runners, runner_out_path)
    if undeut_ground_out_path is not None:
        limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products],
                    undeut_ground_out_path)
    if winner_scores_out_path is not None:
        limit_write(p1.winner_scores, winner_scores_out_path)
    if rtdt_com_cvs_out_path is not None:
        limit_write([p1.rt_com_cv, p1.dt_com_cv], rtdt_com_cvs_out_path)
    if path_plot_out_path is not None:
        undeut_grounds = [p1.undeut_grounds, p1.undeut_ground_dot_products]
        plot_gjr_(winner=p1.winner,
                  undeut_grounds=undeut_grounds,
                  output_path=path_plot_out_path,
                  prefix=name)

    if return_flag:
        out_dict["PathOptimizer"] = p1
        return out_dict

    # Save all ics with all computed attributes to one file
    if not os.path.isdir('resources/ics'):
        os.mkdir('resources/ics')
    limit_write(atc, 'resources/ics/' + name + '.gz.cpickle.zlib')


if __name__ == "__main__":

    # set expected command line arguments
    parser = argparse.ArgumentParser(
        description=
        "Generate a best-estimate HDX-timeseries of IsotopeClusters for a given library protein"
    )
    # inputs
    parser.add_argument("library_info_path", help="path/to/checked_library_info.json")
    parser.add_argument(
        "timepoints_yaml",
        help=
        "path/to/file.yaml containing list of hdx timepoints in integer seconds which are also keys mapping to lists of each timepoint's .mzML file, can pass config/config.yaml - for Snakemake context"
    )
    parser.add_argument(
        "-i",
        "--all_ic_input_paths",
        nargs="*",
        help=
        "structured 2D list of extracted IsotopeCluster objects from each tensor included in the rt_group."
    )
    parser.add_argument(
        "-d",
        "--input_directory_path",
        help=
        "path/to/directory to search for relevant files if assembling filenames automatically, requires --rt_group_name"
    )
    parser.add_argument(
        "-n",
        "--rt_group_name",
        help=
        "rt-group name to use for generating relevant tensor files, requires --input_directory_path"
    )
    parser.add_argument(
        "-g",
        "--old_data_dir",
        help=
        "directory containing Gabe's pickled output files, using this option prints old data on plots"
    )
    # outputs
    parser.add_argument("-o",
                        "--html_plot_out_path",
                        help="path/to/file for .html plot of results")
    parser.add_argument(
        "-w",
        "--winner_out_path",
        help="path/to/file to save winning IsotopeCluster objects")
    parser.add_argument("-r",
                        "--runner_out_path",
                        help="path/to/file to save runner-up IsotopeClusters")
    parser.add_argument(
        "-p",
        "--undeut_ground_out_path",
        help=
        "path/to/file to save selected highest-confidence undeuterated IsotopeClusters"
    )
    parser.add_argument("-s",
                        "--winner_scores_out_path",
                        help="path/to/file to save winning path IC scores")
    parser.add_argument("-c",
                        "--rtdt_com_cvs_out_path",
                        help="path/to/file to save rt/dt error measurement")
    parser.add_argument("-po",
                        "--path_plot_out_path",
                        help="path/to/file to save path plot .pdf")
    parser.add_argument("-mono", "--out_monobody_path", action="store_true",
                        help="generate ic_time_series based on monobody score terms")
    args = parser.parse_args()

    timepoints = yaml.load(open(args.timepoints_yaml, "rb").read(), Loader=yaml.Loader)
    
    # generate explicit inputs and open timpoints .yaml
    if args.all_ic_input_paths is None:
        if args.input_directory_path is not None and args.rt_group_name is not None:
            args.all_ic_input_paths = optimize_paths_inputs(
                args.library_info_path, args.input_directory_path,
                args.rt_group_name, args.timepoints_yaml)
        else:
            parser.print_help()
            sys.exit()


    main(library_info_path=args.library_info_path,
         all_ic_input_paths=args.all_ic_input_paths,
         timepoints=timepoints,
         rt_group_name=args.rt_group_name,
         old_data_dir=args.old_data_dir,
         path_plot_out_path=args.path_plot_out_path,
         html_plot_out_path=args.html_plot_out_path,
         winner_out_path=args.winner_out_path,
         runner_out_path=args.runner_out_path,
         undeut_ground_out_path=args.undeut_ground_out_path,
         winner_scores_out_path=args.winner_scores_out_path,
         rtdt_com_cvs_out_path=args.rtdt_com_cvs_out_path,
         out_monobody_path=args.out_monobody_path)
