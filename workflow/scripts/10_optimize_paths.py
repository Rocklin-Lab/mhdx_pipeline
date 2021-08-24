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
import math
import yaml
import shutil
import argparse
import numpy as np
import pandas as pd
from hdx_limit.core.io import limit_read, limit_write, check_for_create_dirs
from hdx_limit.core.processing import PathOptimizer
from hdx_limit.core.gjr_plot import plot_gjr_


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
    charges = library_info.loc[library_info["name"]==name]['charge'].values
    for key in timepoints["timepoints"]:
        if len(timepoints[key]) > 1:
            for file in timepoints[key]:
                for charge in charges:
                    name_inputs.append(
                        input_directory_path 
                        + name 
                        + "_" 
                        + "charge"
                        + str(charge) 
                        + "_" 
                        + file 
                        +".cpickle.zlib"
                    )  # TODO: This won't work if we go to .mzML/.RAW interoperability.
        else:
            file = timepoints[key][0]
            for charge in charges:
                name_inputs.append(
                    input_directory_path 
                    + name 
                    + "_" 
                    + "charge"
                    + str(charge)
                    + "_"
                    + file 
                    + ".cpickle.zlib")

    return name_inputs

def gen_correlate_matrix(list_of_arrays):
    """Takes a list of 1D signal arrays and generates a correlation matrix for those signals. 
    Args:
        list_of_arrays (list of number arrays): A flat iterable of signals (arrays) to be compared.

    Returns:
        corr_matrix (numpy array): A correlation matrix containing correlation coefficients for each pair of signal arrays.

    """
    corr_matrix = np.zeros((len(list_of_arrays), len(list_of_arrays)))
    for ind1, arr1 in enumerate(list_of_arrays):
        for ind2, arr2 in enumerate(list_of_arrays):
            corr_matrix[ind1, ind2] = max(np.correlate(arr1 / np.linalg.norm(arr1), arr2 / np.linalg.norm(arr2)))

    return corr_matrix


def main(library_info_path,
         all_ic_input_paths,
         timepoints,
         monobody_return_flag=False,
         multibody_return_flag=False,
         rt_group_name=None,
         old_data_dir=None,
         all_timepoint_clusters_out_path=None,
         prefiltered_ics_out_path=None,
         mono_path_plot_out_path=None,
         mono_html_plot_out_path=None,
         mono_winner_out_path=None,
         mono_runner_out_path=None,
         mono_undeut_ground_out_path=None,
         mono_winner_scores_out_path=None,
         mono_rtdt_com_cvs_out_path=None,
         multi_path_plot_out_path=None,
         multi_html_plot_out_path=None,
         multi_winner_out_path=None,
         multi_runner_out_path=None,
         multi_undeut_ground_out_path=None,
         multi_winner_scores_out_path=None,
         multi_rtdt_com_cvs_out_path=None):
    """Uses PathOptimzier class to generate best-estimate hdx-timeseries of IsotopeClusters for a given library protein.

    Args:
        library_info_path (str): path/to/checked_library_info.json
        all_ic_input_paths (list of strings): list of paths/to/files.cpickle.zlib for all lists of IsotopeClusters from generate_tensor_ics.py
        timepoints (dict): dictionary with 'timepoints' key containing list of hdx timepoints in integer seconds, which are keys mapping to lists of each timepoint's replicate .mzML filenames
        monobody_return_flag: option to return monobody output in python, for notebook context, can be combined with multibody_return_flag.
        multibody_return_flag: option to return multibody output in python, for notebook context, can be combined with monobody_return_flag.
        rt_group_name (str): library_info['name'] value
        old_data_dir (str): path/to/dir to provide comparison to GJR formatted results
        all_timepoint_clusters_out_path (str): path/to/file to output all clusters collected for PathOptimizer as a nested list.
        prefiltered_ics_out_path (str): path/to/file to output ICs selected from prefiltering as a nested list.
        mono_path_plot_out_path (str): path/to/dir for pdf plot of monobody scoring results
        mono_html_plot_out_path (str): path/to/file.html for interactive bokeh plot of monobody scoring results
        mono_winner_out_path (str): path/to/file for winning path from monobody scoring
        mono_runner_out_path (str): path/to/file for top n_runners paths from monobody scoring
        mono_undeut_ground_out_path (str): path/to/file for undeuterated ground-truth IsotopeClusters from monobody scoring
        mono_winner_scores_out_path (str): path/to/file for monobody scoring winning path score values 
        mono_rtdt_com_cvs_out_path (str): path/to/file for rt and dt correlation values from monobody scoring
        multi_path_plot_out_path (str): path/to/dir for pdf plot of multibody scoring results
        multi_html_plot_out_path (str): path/to/file.html for interactive bokeh plot of multibody results
        multi_winner_out_path (str): path/to/file for winning path from multibody scoring
        multi_runner_out_path (str): path/to/file for top n_runners paths from multibody scoring
        multi_undeut_ground_out_path (str): path/to/file for undeuterated ground-truth IsotopeClusters from multibody scoring
        multi_winner_scores_out_path (str): path/to/file for multibody winning path score values 
        multi_rtdt_com_cvs_out_path (str): path/to/file for rt and dt correlation values from multibody scoring

    Returns:
        out_dict (dict): dictionary containing 'path_optimizer' key, and corresponding PathOptimizer object 

    """
    monobody_path_arguments = [
                                    mono_path_plot_out_path,
                                    mono_html_plot_out_path,
                                    mono_winner_out_path, 
                                    mono_runner_out_path, 
                                    mono_undeut_ground_out_path, 
                                    mono_winner_scores_out_path, 
                                    mono_rtdt_com_cvs_out_path,
                                ]

    multibody_path_arguments = [
                                    multi_path_plot_out_path,
                                    multi_html_plot_out_path,
                                    multi_winner_out_path, 
                                    multi_runner_out_path, 
                                    multi_undeut_ground_out_path, 
                                    multi_winner_scores_out_path, 
                                    multi_rtdt_com_cvs_out_path,
                                ]

    out_dict = {}
    library_info = pd.read_json(library_info_path)

    if rt_group_name is None:
        name = all_ic_input_paths[0].split("/")[-2]
    else:
        name = rt_group_name

    # Order files, pooling all replicates and charges by timepoint.
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

    # Generate nearest_neighbor_correlation for each ic.
    # TODO: Maybe atc should be saved as a gz.cpickle.zlib file?
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

    # Starting condition output arguments.
    if all_timepoint_clusters_out_path is not None:
        limit_write(atc, all_timepoint_clusters_out_path)
    if prefiltered_ics_out_path is not None:
        limit_write(p1.prefiltered_ics, prefiltered_ics_out_path)
    
    # Checks if arguments require monobody scoring run.
    if (any(arg is not None for arg in monobody_path_arguments)) or (monobody_return_flag is not False):

        p1.optimize_paths_mono()

        if monobody_return_flag is not False:
            out_dict["monobody_winner"] = p1.winner
            out_dict["monobody_runners"] = p1.runners
            out_dict["monobody_undeut_grounds"] = [p1.undeut_grounds, p1.undeut_ground_dot_products]
            out_dict["monobody_winner_scores"] = p1.winner_scores
            out_dict["monobody_rt_dt_com_cvs"] = [p1.rt_com_cv, p1.dt_com_cv]

        if mono_path_plot_out_path is not None:
            undeut_grounds = [p1.undeut_grounds, p1.undeut_ground_dot_products]
            plot_gjr_(winner=p1.winner,
                      undeut_grounds=undeut_grounds,
                      output_path=mono_path_plot_out_path,
                      prefix=name)
        if mono_html_plot_out_path is not None:
             p1.bokeh_plot(mono_html_plot_out_path)
        if mono_winner_out_path is not None:
            limit_write(p1.winner, mono_winner_out_path)
        if mono_runner_out_path is not None:
            limit_write(p1.runners, mono_runner_out_path)
        if mono_undeut_ground_out_path is not None:
            limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products],
                        mono_undeut_ground_out_path)
        if mono_winner_scores_out_path is not None:
            limit_write(p1.winner_scores, mono_winner_scores_out_path)
        if mono_rtdt_com_cvs_out_path is not None:
            limit_write([p1.rt_com_cv, p1.dt_com_cv], mono_rtdt_com_cvs_out_path)

    # Checks if arguments require multibody scoring run.
    if (any(arg is not None for arg in multibody_path_arguments)) or (multibody_return_flag is not False):
       
        p1.optimize_paths_multi()

        if monobody_return_flag is not False:
            out_dict["monobody_winner"] = p1.winner
            out_dict["monobody_runners"] = p1.runners
            out_dict["monobody_undeut_grounds"] = [p1.undeut_grounds, p1.undeut_ground_dot_products]
            out_dict["monobody_winner_scores"] = p1.winner_scores
            out_dict["monobody_rt_dt_com_cvs"] = [p1.rt_com_cv, p1.dt_com_cv]

        if multi_path_plot_out_path is not None:
            undeut_grounds = [p1.undeut_grounds, p1.undeut_ground_dot_products]
            plot_gjr_(winner=p1.winner,
                      undeut_grounds=undeut_grounds,
                      output_path=multi_path_plot_out_path,
                      prefix=name)
        if multi_html_plot_out_path is not None:
             p1.bokeh_plot(multi_html_plot_out_path)
        if multi_winner_out_path is not None:
            limit_write(p1.winner, multi_winner_out_path)
        if multi_runner_out_path is not None:
            limit_write(p1.runners, multi_runner_out_path)
        if multi_undeut_ground_out_path is not None:
            limit_write([p1.undeut_grounds, p1.undeut_ground_dot_products],
                        multi_undeut_ground_out_path)
        if multi_winner_scores_out_path is not None:
            limit_write(p1.winner_scores, multi_winner_scores_out_path)
        if multi_rtdt_com_cvs_out_path is not None:
            limit_write([p1.rt_com_cv, p1.dt_com_cv], multi_rtdt_com_cvs_out_path)


if __name__ == "__main__":

    # Checks for Snakemake context and parses arguments.
    if "snakemake" in globals():

        library_info_path = snakemake.input[0]
        timepoints = yaml.load(open(snakemake.input[1], "rb").read(), Loader=yaml.Loader)
        all_ic_input_paths = snakemake.input[2:]
        old_data_dir = None
        rt_group_name = snakemake.params.rt_group_name
        all_timepoint_clusters_out_path = snakemake.output[0]
        prefiltered_ics_out_path = snakemake.output[1]
        mono_path_plot_out_path = snakemake.output[2]
        mono_html_plot_out_path = None
        mono_winner_out_path = snakemake.output[3]
        mono_runner_out_path = snakemake.output[4]
        mono_undeut_ground_out_path = snakemake.output[5]
        mono_winner_scores_out_path = snakemake.output[6]
        mono_rtdt_com_cvs_out_path = snakemake.output[7]
        multi_path_plot_out_path = snakemake.output[8]
        multi_html_plot_out_path = None
        multi_winner_out_path = snakemake.output[9]
        multi_runner_out_path = snakemake.output[10]
        multi_undeut_ground_out_path = snakemake.output[11]
        multi_winner_scores_out_path = snakemake.output[12]
        multi_rtdt_com_cvs_out_path = snakemake.output[13]

        main(library_info_path=library_info_path,
             all_ic_input_paths=all_ic_input_paths,
             timepoints=timepoints,
             old_data_dir = old_data_dir,
             rt_group_name=rt_group_name,
             all_timepoint_clusters_out_path=all_timepoint_clusters_out_path,
             prefiltered_ics_out_path=prefiltered_ics_out_path,
             mono_path_plot_out_path=mono_path_plot_out_path,
             mono_html_plot_out_path=mono_html_plot_out_path,
             mono_winner_out_path=mono_winner_out_path,
             mono_runner_out_path=mono_runner_out_path,
             mono_undeut_ground_out_path=mono_undeut_ground_out_path,
             mono_winner_scores_out_path=mono_winner_scores_out_path,
             mono_rtdt_com_cvs_out_path=mono_rtdt_com_cvs_out_path,
             multi_path_plot_out_path=multi_path_plot_out_path,
             multi_html_plot_out_path=multi_html_plot_out_path,
             multi_winner_out_path=multi_winner_out_path,
             multi_runner_out_path=multi_runner_out_path,
             multi_undeut_ground_out_path=multi_undeut_ground_out_path,
             multi_winner_scores_out_path=multi_winner_scores_out_path,
             multi_rtdt_com_cvs_out_path=multi_rtdt_com_cvs_out_path)

    else:

        # Configure command line argument parser.
        parser = argparse.ArgumentParser(
            description=
            "Generate a best-estimate HDX-timeseries of IsotopeClusters for a given library protein"
        )
        # Inputs
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
            "directory containing Gabe's pickled output files, using this option prints old data on html plots"
        )
        
        # Starting condition output arguments.
        parser.add_argument("--all_timepoint_clusters_out_path",
                            help="path/to/file to output all input IsotopeClusters as a nested list")
        parser.add_argument("--prefiltered_ics_out_path",
                            help="path/to/file to output ICs passing prefiltering as a nested list")

        # Monobody scoring outputs.
        parser.add_argument("--mono_html_plot_out_path",
                            help="path/to/file for .html plot of results")
        parser.add_argument(
            "--mono_winner_out_path",
            help="path/to/file to save winning IsotopeCluster objects")
        parser.add_argument("--mono_runner_out_path",
                            help="path/to/file to save runner-up IsotopeClusters")
        parser.add_argument(
            "--mono_undeut_ground_out_path",
            help=
            "path/to/file to save selected highest-confidence undeuterated IsotopeClusters"
        )
        parser.add_argument("--mono_winner_scores_out_path",
                            help="path/to/file to save winning path IC scores")
        parser.add_argument("--mono_rtdt_com_cvs_out_path",
                            help="path/to/file to save rt/dt error measurement")
        parser.add_argument("--mono_path_plot_out_path",
                            help="path/to/file to save path plot .pdf")
        
        # Multibody scoring outputs.
        parser.add_argument("--multi_html_plot_out_path",
                            help="path/to/file for .html plot of results")
        parser.add_argument(
            "--multi_winner_out_path",
            help="path/to/file to save winning IsotopeCluster objects")
        parser.add_argument("--multi_runner_out_path",
                            help="path/to/file to save runner-up IsotopeClusters")
        parser.add_argument(
            "--multi_undeut_ground_out_path",
            help=
            "path/to/file to save selected highest-confidence undeuterated IsotopeClusters"
        )
        parser.add_argument("--multi_winner_scores_out_path",
                            help="path/to/file to save winning path IC scores")
        parser.add_argument("--multi_rtdt_com_cvs_out_path",
                            help="path/to/file to save rt/dt error measurement")
        parser.add_argument("--multi_path_plot_out_path",
                            help="path/to/file to save path plot .pdf")
        
        args = parser.parse_args()

        # Opens timepoints .yaml and generates explicit inputs.
        timepoints = yaml.load(open(args.timepoints_yaml, "rb").read(), Loader=yaml.Loader)
        if args.all_ic_input_paths is None:
            if args.input_directory_path is not None and args.rt_group_name is not None:
                args.all_ic_input_paths = optimize_paths_inputs(
                    args.library_info_path, args.input_directory_path,
                    args.rt_group_name, timepoints)
            else:
                parser.print_help()
                sys.exit()

        # Lists of arguments grouped by required subdirectory.
        output_paths = [    
                            args.all_timepoint_clusters_out_path,
                            args.prefiltered_ics_out_path,
                            args.mono_winner_out_path, 
                            args.mono_runner_out_path, 
                            args.mono_undeut_ground_out_path, 
                            args.mono_winner_scores_out_path, 
                            args.mono_rtdt_com_cvs_out_path,
                            args.mono_path_plot_out_path,
                            args.mono_html_plot_out_path,
                            args.multi_winner_out_path, 
                            args.multi_runner_out_path, 
                            args.multi_undeut_ground_out_path, 
                            args.multi_winner_scores_out_path, 
                            args.multi_rtdt_com_cvs_out_path,
                            args.multi_path_plot_out_path,
                            args.multi_html_plot_out_path,
                        ]

        # Checks for arguments that require additional directories and creates them if they don't exist.
        check_for_create_dirs(output_paths)

        main(library_info_path=args.library_info_path,
             all_ic_input_paths=args.all_ic_input_paths,
             timepoints=timepoints,
             rt_group_name=args.rt_group_name,
             old_data_dir=args.old_data_dir,
             all_timepoint_clusters_out_path=args.all_timepoint_clusters_out_path,
             prefiltered_ics_out_path=args.prefiltered_ics_out_path,
             mono_path_plot_out_path=args.mono_path_plot_out_path,
             mono_html_plot_out_path=args.mono_html_plot_out_path,
             mono_winner_out_path=args.mono_winner_out_path,
             mono_runner_out_path=args.mono_runner_out_path,
             mono_undeut_ground_out_path=args.mono_undeut_ground_out_path,
             mono_winner_scores_out_path=args.mono_winner_scores_out_path,
             mono_rtdt_com_cvs_out_path=args.mono_rtdt_com_cvs_out_path,
             multi_path_plot_out_path=args.multi_path_plot_out_path,
             multi_html_plot_out_path=args.multi_html_plot_out_path,
             multi_winner_out_path=args.multi_winner_out_path,
             multi_runner_out_path=args.multi_runner_out_path,
             multi_undeut_ground_out_path=args.multi_undeut_ground_out_path,
             multi_winner_scores_out_path=args.multi_winner_scores_out_path,
             multi_rtdt_com_cvs_out_path=args.multi_rtdt_com_cvs_out_path)
