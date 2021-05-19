import os
import sys
import yaml
import psutil
import numpy as np
import argparse
import pandas as pd

sys.path.append(os.getcwd() + "/workflow/scripts/")
from HDX_LIMIT.processing import TensorGenerator, generate_tensor_factors
from HDX_LIMIT.io import limit_write


def main(library_info_path,
         tensor_input_path,
         timepoints_dict,
         isotope_clusters_out_path=None,
         factor_out_path=None,
         factor_plot_output_path=None,
         return_flag=False,
         gauss_params=(3, 1),
         n_factors=15,
         filter_factors=False,
         factor_rt_r2_cutoff=0.91,
         factor_dt_r2_cutoff=0.91,
         ic_peak_prominence=0.15,
         ic_peak_width=3,
         ic_rel_height_filter=True,
         ic_rel_height_filter_baseline=0.15,
         ic_rel_height_threshold=0.10):
    """Performs nonnegative tensor factorization to deconvolute input tensor, identifies IsotopeCluster objects, 
    and optionally returns or writes output list of IsotopeClusters.

    Args:
        library_info_path (str): path/to/library_info.csv
        tensor_input_path (str): path/to/tensor.cpickle.zlib
        timepoints_dict (dict): dictionary with 'timepoints' key containing list of hdx timepoints in integer seconds, which are keys mapping to lists of each timepoint's replicate .mzML filenames 
        isotope_clusters_out_path (str): path/to/file for main output - list of IsotopeClusters objects
        return_flag (bool): option to return output in python, for notebook context
        gauss_params (tuple of ints/floats): Gaussian smoothing parameters in LC-RT and IMS-DT dimensions, (rt_sigma, dt_sigma)

    Returns:
        out_dict (dict): dictionary containing TensorGenerator object
    
    """
    out_dict = {}

    # open library_info
    library_info = pd.read_csv(library_info_path)
    my_idx = int(tensor_input_path.split("/")[-1].split("_")[0])
    # my_centers = library_info.iloc[my_idx]["mz_centers"].values
    ## temporary fix for getting mz center values in an array
    my_centers = library_info.iloc[my_idx]["mz_centers"]
    centers = my_centers.split()
    centers = [x.replace("[", "") for x in centers]
    centers = [x.replace("]", "") for x in centers]
    centers = np.array([float(x) for x in centers])

    # find timepoint of passed filename by config comparison
    for tp in timepoints_dict["timepoints"]:
        for fn in timepoints_dict[tp]:
            if fn in tensor_input_path:
                my_tp = tp

    data_tensor = generate_tensor_factors(tensor_fpath=tensor_input_path,
                                          library_info_df=library_info,
                                          timepoint_index=my_tp,
                                          gauss_params=gauss_params,
                                          n_factors=n_factors,
                                          mz_centers=centers,
                                          factor_output_fpath=factor_out_path,
                                          factor_plot_output_path=factor_plot_output_path,
                                          timepoint_label=None,
                                          filter_factors=filter_factors,
                                          factor_rt_r2_cutoff=factor_rt_r2_cutoff,
                                          factor_dt_r2_cutoff=factor_dt_r2_cutoff)

    all_ics = []
    for factor in data_tensor.DataTensor.factors:

        # generate isotope cluster class
        factor.find_isotope_clusters(prominence=ic_peak_prominence,
                                     width_val=ic_peak_width,
                                     rel_height_filter=ic_rel_height_filter,
                                     baseline_threshold=ic_rel_height_filter_baseline,
                                     rel_height_threshold=ic_rel_height_threshold)

        for ic in factor.isotope_clusters:
            all_ics.append(ic)

    if isotope_clusters_out_path is not None:
        limit_write(all_ics, isotope_clusters_out_path)

    if return_flag:
        out_dict["TensorGenerator"] = data_tensor
        return out_dict


if __name__ == "__main__":

    # set expected command line arguments
    parser = argparse.ArgumentParser(
        description=
        "Accepts tensor as input, factorizes and saves IsotopeClusters from resulting Factors"
    )
    parser.add_argument("library_info_path", help="path/to/library_info.csv")
    parser.add_argument(
        "tensor_input_path",
        help="path/to/file.cpickle.zlib for tensor to factorize")
    parser.add_argument(
        "timepoints_yaml",
        help=
        "path/to/file.yaml containing list of hdx timepoints in integer seconds which are also keys mapping to lists of each timepoint's .mzML file, can pass config/config.yaml - for Snakemake context"
    )
    parser.add_argument(
        "-o",
        "--isotope_clusters_out_path",
        help="path/to/output.cpickle.zlib, list of IsotopeClusters")
    parser.add_argument(
        "-of",
        "--factor_data_out_path",
        help="path/to/output.cpickle.zlib.factor, FactorData")
    parser.add_argument(
        "-po",
        "--factor_plot_out_path",
        help="path/to/output.cpickle.zlib.factor.pdf, FactorData Plot output .pdf")
    parser.add_argument(
        "-r",
        "--return_flag",
        type=bool,
        default=False,
        help="option to return output dictionary in python, for notebook context"
    )
    parser.add_argument(
        "-g",
        "--gauss_params",
        type=tuple,
        default=(3, 1),
        help="determines intensity of gaussian smoothing in rt and dt dimensions"
    )
    parser.add_argument(
        "-n",
        "--n_factors",
        type=int,
        default=15,
        help="maximum number of factors for factorization of the data tensor"
    )
    args = parser.parse_args()

    # open timepoints .yaml into dict for main()
    config_dict = yaml.load(open(args.timepoints_yaml, 'rb'), Loader=yaml.Loader)

    filter_factors = config_dict["filter_factor"]
    factor_rt_r2_cutoff = config_dict["factor_rt_r2_cutoff"]
    factor_dt_r2_cutoff = config_dict["factor_dt_r2_cutoff"]

    ic_peak_prom = config_dict["ic_peak_prominence"]
    ic_peak_width = config_dict["ic_peak_width"]
    ic_rel_ht_filter = config_dict["ic_rel_height_filter"]
    ic_rel_ht_baseline = config_dict["ic_rel_height_filter_baseline"]
    ic_rel_ht_threshold = config_dict["ic_rel_height_threshold"]



    main(library_info_path=args.library_info_path,
         tensor_input_path=args.tensor_input_path,
         timepoints_dict=config_dict,
         isotope_clusters_out_path=args.isotope_clusters_out_path,
         factor_out_path=args.factor_data_out_path,
         factor_plot_output_path=args.factor_plot_out_path,
         return_flag=args.return_flag,
         gauss_params=args.gauss_params,
         filter_factors=filter_factors,
         factor_rt_r2_cutoff=factor_rt_r2_cutoff,
         factor_dt_r2_cutoff=factor_dt_r2_cutoff,
         ic_peak_prominence=ic_peak_prom,
         ic_peak_width=ic_peak_width,
         ic_rel_height_filter=ic_rel_ht_filter,
         ic_rel_height_filter_baseline=ic_rel_ht_baseline,
         ic_rel_height_threshold=ic_rel_ht_threshold)
