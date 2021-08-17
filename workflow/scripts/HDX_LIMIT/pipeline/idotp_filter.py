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
import sys
import glob
import argparse
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


def main(library_info_path,
         all_idotp_csv_inputs,
         indices_out_path=None,
         library_info_out_path=None,
         plot_out_path=None,
         return_flag=False,
         idotp_cutoff=0.95):
    """Reads all library_info index idotp_check.csv files and returns or saves a list of indices with idotp >= idotp_cutoff.

    Args:
        library_info_path (str): path/to/library_info.json
        all_idotp_csv_inputs (list of strings): list of all input IsotopeCluster-list filepaths
        indices_out_path (str): path/to/filter_passing_indices.csv
        library_info_out_path (str): path/to/checked_library_info.json
        plot_out_path (str): path/to/file.png for idotp_distribution plot
        return_flag (bool): option to return a dictionary of outputs in python context
        idotp_cutoff (float): inclusive lower-bound on idotp [0,1] to be considered for evaluation, default=0.95

    Returns:
        out_dict (dict) = dictionary containing "filter_passing_indices"

    """
    library_info = pd.read_json(library_info_path)
    sorted_inputs = sorted(all_idotp_csv_inputs, key=lambda fn: int(fn.split("_")[-3]))
    print("Length of inputs: "+str(len(sorted_inputs)))

    out_dict = {}
    filter_passing_indices = []
    idotps = []
    mz_centers = []
    theor_mz_dists = []
    integrated_mz_width_list = []

    for fn in sorted_inputs:
        prot_name = fn.split("/")[-2] # Name from protein directory.
        prot_charge = int([item[6:] for item in fn.split("/")[-1].split("_") if "charge" in item][0]) # Finds by keyword and strip text.
        lib_idx = library_info.loc[(library_info["name"]==prot_name) & (library_info["charge"]==prot_charge)].index
        idpc = pd.read_json(fn)
        idotps.append(idpc["idotp"].values[0])
        mz_centers.append(idpc["mz_centers"][0]) # Account for nested list structure
        theor_mz_dists.append(idpc["theor_mz_dist"][0])
        integrated_mz_width_list.append(idpc["integrated_mz_width"].values[0])
        if idpc["idotp"].values[0] >= idotp_cutoff:
            filter_passing_indices.append(lib_idx)

    # Set values in library_info and write out
    library_info["idotp"] = idotps
    library_info["mz_centers"] = mz_centers
    library_info["theor_mz_dist"] = theor_mz_dists
    library_info["integrated_mz_width"] = integrated_mz_width_list

    if library_info_out_path is not None:
        library_info.to_json(library_info_out_path)

    # re-order indices
    filter_passing_indices = sorted(filter_passing_indices)
    # add passing indices to output dict

    out_dict["filter_passing_indices"] = filter_passing_indices
    out_dict["mz_centers"] = mz_centers
    out_dict["theor_mz_dist"] = theor_mz_dists

    out_df = pd.DataFrame.from_dict({"index": filter_passing_indices})

    if plot_out_path is not None:
        sns.displot(idotps)
        plt.axvline(idotp_cutoff, 0, 1)
        plt.savefig(plot_out_path)

    if indices_out_path is not None:
        out_df.to_csv(indices_out_path)

    if return_flag:
        return out_dict


if __name__ == "__main__":
    # If the snakemake global object is present, save expected arguments from snakemake to be passed to main().
    if "snakemake" in globals():
        library_info_path = snakemake.input.pop(0)
        all_idotp_csv_inputs = snakemake.input

        indices_out_path = snakemake.output[0]
        library_info_out_path =  snakemake.output[1]
        plot_out_path = snakemake.output[2]


        main(library_info_path=library_info_path,
             all_idotp_csv_inputs=all_idotp_csv_inputs,
             indices_out_path=indices_out_path,
             library_info_out_path=library_info_out_path,
             plot_out_path=plot_out_path)
    else:
        # CLI context, set expected arguments with argparse module.
        parser = argparse.ArgumentParser(
            description=
            "Reads all rt-group idotp csvs and returns or saves a list of indices with idotp >= idotp_cutoff."
        )
        parser.add_argument("library_info_path", help="path/to/library_info.json")
        parser.add_argument("-i",
                            "--all_idotp_csv_inputs",
                            help="list of all idotp check .csv outputs to be read")
        parser.add_argument("-d",
                            "--input_dir_path",
                            help="path/to/dir/ containing idotp_check.csv files")
        parser.add_argument("-o",
                            "--indices_out_path",
                            help="path/to/filter_passing_indices.csv")
        parser.add_argument("-l",
                            "--library_info_out_path",
                            help="path/to/checked_library_info.json")
        parser.add_argument("--p",
                            "--plot_out_path",
                            help="path/to/idotp_distribution.png")
        parser.add_argument(
            "-c",
            "--idotp_cutoff",
            type=float,
            default=0.95,
            help=
            "lower limit on dot-product between theoretical integrated m/z of POI and int. m/z of observed signal in question. Float in range [0,1], default 0.95 "
        )
        args = parser.parse_args()

        if args.all_idotp_csv_inputs is None and args.input_dir_path is None:
            parser.print_help()
            sys.exit()

        if args.all_idotp_csv_inputs is None and args.input_dir_path is not None:
            args.all_idotp_csv_inputs = sorted(
                list(glob.glob(args.input_dir_path + "*idotp_check.csv")))

        all_idotp_csv_inputs = args.all_idotp_csv_inputs.split(' ')

        main(library_info_path=args.library_info_path,
             all_idotp_csv_inputs=args.all_idotp_csv_inputs,
             indices_out_path=args.indices_out_path,
             library_info_out_path = args.library_info_out_path,
             plot_out_path=args.plot_out_path,
             idotp_cutoff=args.idotp_cutoff)
