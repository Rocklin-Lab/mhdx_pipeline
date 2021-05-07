import sys
import glob
import argparse
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


def main(all_idotp_csv_inputs,
         out_path=None,
         plot_out_path=None,
         return_flag=False,
         idotp_cutoff=0.95):
    """Reads all library_info index idotp_check.csv files and returns or saves a list of indices with idotp >= idotp_cutoff.

    Args:
        all_idotp_csv_inputs (list of strings): list of all input IsotopeCluster-list filepaths
        out_path (str): path/to/file for main output.csv
        plot_out_path (str): path/to/file for idotp_distribution plot
        return_flag (bool): option to return main output in python, for notebook context
        idotp_cutoff (float): inclusive lower-bound on idotp [0,1] to be considered for evaluation, default=0.95

    Returns:
        out_dict (dict) = dictionary containing "filter_passing_indices"

    """
    out_dict = {}

    filter_passing_indices = []
    idotps = []
    for fn in all_idotp_csv_inputs:
        lib_idx = int(fn.split("/")[-1].split("_")[0])
        idpc = pd.read_csv(fn)
        idotps.append(idpc["idotp"].values[0])
        if idpc["idotp"].values[0] >= idotp_cutoff:
            filter_passing_indices.append(lib_idx)

    # re-order indices
    filter_passing_indices = sorted(filter_passing_indices)
    # add passing indices to output dict
    out_dict["index"] = filter_passing_indices
    # make df output option
    out_df = pd.DataFrame.from_dict(out_dict)

    if plot_out_path is not None:
        sns.displot(idotps)
        plt.savefig(plot_out_path)

    if out_path is not None:
        out_df.to_csv(out_path)

    if return_flag:
        return out_dict


if __name__ == "__main__":

    # set expected command line arguments
    parser = argparse.ArgumentParser(
        description=
        "Reads all rt-group idotp csvs and returns or saves a list of indices with idotp >= idotp_cutoff."
    )
    parser.add_argument("-i",
                        "--all_idotp_csv_inputs",
                        help="list of all idotp check .csv outputs to be read")
    parser.add_argument("-d",
                        "--input_dir_path",
                        help="path/to/dir/ containing idotp_check.csv files")
    parser.add_argument("-o",
                        "--out_path",
                        help="path/to/filter_passing_indices.csv")
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

    main(args.all_idotp_csv_inputs,
         out_path=args.out_path,
         idotp_cutoff=args.idotp_cutoff)
