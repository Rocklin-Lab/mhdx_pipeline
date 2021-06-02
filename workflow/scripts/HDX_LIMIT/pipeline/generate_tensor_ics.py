import os
import sys
import yaml
import psutil
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

sys.path.append(os.getcwd() + "/workflow/scripts/")
from HDX_LIMIT.processing import TensorGenerator, generate_tensor_factors
from HDX_LIMIT.io import limit_write



def convert_factor_mz_data_to_cont(mz_data, mz_labels, bins_per_isotope_peak):
    padded_mz_labels = []
    padded_factor_mz = []
    mz_label_spacing = mz_labels[1] - mz_labels[0]

    for i in range(0, len(mz_labels), bins_per_isotope_peak):
        padded_mz_labels.append(mz_labels[i] - mz_label_spacing)
        padded_mz_labels += list(mz_labels[i:i + bins_per_isotope_peak])
        padded_mz_labels.append(padded_mz_labels[-1] + mz_label_spacing)

        padded_factor_mz.append(0)
        padded_factor_mz += list(mz_data[i:i + bins_per_isotope_peak])
        padded_factor_mz.append(0)

    return padded_mz_labels, padded_factor_mz


def plot_mz_data(fig, gs, row_num, col_num, mz_label, mz_data, plot_label):

    mz_sum = np.round(np.sum(mz_data), 2)

    ax = fig.add_subplot(gs[row_num, col_num])
    plt.plot(mz_label, mz_data, linewidth=0.5, label=plot_label)
    plt.text(0, 0.1, s=str(mz_sum), fontsize=10, transform=ax.transAxes)
    ax.set_yticks([])
    ax.tick_params(length=3, pad=3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.legend(loc='best', fontsize='small')


def plot_ics(list_of_ics_from_a_factor):


    num_of_mz_plots = len(list_of_ics_from_a_factor) + 1

    num_columns = 3
    num_rows = 0

    for num in range(num_of_mz_plots):
        if num % num_columns == 0:
            num_rows += 1


    pad_factor_mz_label, pad_factor_mz = convert_factor_mz_data_to_cont(
        mz_data=list_of_ics_from_a_factor[0].factor_mz_data,
        mz_labels=list_of_ics_from_a_factor[0].mz_labels,
        bins_per_isotope_peak=list_of_ics_from_a_factor[0].bins_per_isotope_peak)

    fig = plt.figure(figsize=(15, num_rows * 1.6))
    gs = gridspec.GridSpec(ncols=num_columns, nrows=num_rows, figure=fig)

    n_rows = 0
    n_cols = 0

    # plot factor mz data
    plot_mz_data(fig=fig,
                 gs=gs,
                 row_num=n_rows,
                 col_num=n_cols,
                 mz_label=pad_factor_mz_label,
                 mz_data=pad_factor_mz,
                 plot_label='factor_mz')


    # plot ics from the factor
    n_cols = 1

    for num, ic in enumerate(list_of_ics_from_a_factor):

        ic_mz_label, ic_mz_ = convert_factor_mz_data_to_cont(mz_data=ic.cluster_mz_data,
                                                             mz_labels=ic.mz_labels,
                                                             bins_per_isotope_peak=ic.bins_per_isotope_peak)

        plot_mz_data(fig=fig,
                     gs=gs,
                     row_num=n_rows,
                     col_num=n_cols,
                     mz_label=ic_mz_label,
                     mz_data=ic_mz_,
                     plot_label=str(num))

        if (n_cols+1) % num_columns == 0:
            n_rows += 1
            n_cols = 0
        else:
            n_cols += 1

    title = 'Factor num %s' % list_of_ics_from_a_factor[0].factor_idx

    plt.suptitle(title)
    plt.tight_layout()


def plot_ics_from_ic_list(list_of_ics, output_path):

    ### sort ics with factor number
    factor_indices = [x.factor_idx for x in list_of_ics]
    unique_factor_idx = np.unique(factor_indices)
    sorted_list_ics = [[] for _ in range(len(unique_factor_idx))]


    for ind, factor_idx in enumerate(unique_factor_idx):
        for ics in list_of_ics:
            if ics.factor_idx == factor_idx:
                sorted_list_ics[ind].append(ics)

    ### plotting ics
    with PdfPages(output_path) as pdf:
        for ind, ic_list in enumerate(sorted_list_ics):
            plot_ics(ic_list)

            pdf.savefig()
            plt.close()


def main(library_info_path,
         tensor_input_path,
         timepoints_dict,
         normalization_factors,
         isotope_clusters_out_path=None,
         factor_out_path=None,
         factor_plot_output_path=None,
         ic_plot_output_path=None,
         return_flag=False,
         gauss_params=(3, 1),
         n_factors=15,
         filter_factors=False,
         factor_rt_r2_cutoff=0.91,
         factor_dt_r2_cutoff=0.91,
         ic_peak_prominence=0.10,
         auto_ic_peak_width=True,
         ic_peak_width=2,
         ic_rel_height_filter=True,
         ic_rel_height_filter_baseline=0.10,
         ic_rel_height_threshold=0.10):
    """Performs nonnegative tensor factorization to deconvolute input tensor, identifies IsotopeCluster objects, 
    and optionally returns or writes output list of IsotopeClusters.

    Args:
        library_info_path (str): path/to/library_info.json
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
    library_info = pd.read_json(library_info_path)
    my_idx = int(tensor_input_path.split("/")[-1].split("_")[0])
    # my_centers = library_info.iloc[my_idx]["mz_centers"].values
    ## temporary fix for getting mz center values in an array
    my_centers = library_info.iloc[my_idx]["mz_centers"]
    centers = np.array(my_centers)

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
                                          normalization_factors=normalization_factors,
                                          factor_output_fpath=factor_out_path,
                                          factor_plot_output_path=factor_plot_output_path,
                                          timepoint_label=None,
                                          filter_factors=filter_factors,
                                          factor_rt_r2_cutoff=factor_rt_r2_cutoff,
                                          factor_dt_r2_cutoff=factor_dt_r2_cutoff)

    all_ics = []

    ic_peak_width_auto = 0.8 * library_info['integrated_mz_width'].values[my_idx]

    for factor in data_tensor.DataTensor.factors:

        # generate isotope cluster class

        if auto_ic_peak_width:
            factor.find_isotope_clusters(prominence=ic_peak_prominence,
                                     width_val=ic_peak_width_auto,
                                     rel_height_filter=ic_rel_height_filter,
                                     baseline_threshold=ic_rel_height_filter_baseline,
                                     rel_height_threshold=ic_rel_height_threshold)
        else:
            factor.find_isotope_clusters(prominence=ic_peak_prominence,
                                         width_val=ic_peak_width,
                                         rel_height_filter=ic_rel_height_filter,
                                         baseline_threshold=ic_rel_height_filter_baseline,
                                         rel_height_threshold=ic_rel_height_threshold)

        for ic in factor.isotope_clusters:
            all_ics.append(ic)

    if isotope_clusters_out_path is not None:
        limit_write(all_ics, isotope_clusters_out_path)

    if ic_plot_output_path is not None:
        plot_ics_from_ic_list(all_ics, ic_plot_output_path)

    if return_flag:
        out_dict["TensorGenerator"] = data_tensor
        return out_dict


if __name__ == "__main__":

    # set expected command line arguments
    parser = argparse.ArgumentParser(
        description=
        "Accepts tensor as input, factorizes and saves IsotopeClusters from resulting Factors"
    )
    parser.add_argument("library_info_path", help="path/to/library_info.json")
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
        "-po_ic",
        "--ic_plot_out_path",
        help="path/to/output.cpickle.zlib.ic.pdf, ICData Plot output .pdf")
    parser.add_argument(
        "-p",
        "--normalization_factors_path",
        help="path/to/normalization_factors.csv")
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
    auto_ic_peak_width = config_dict["auto_ic_peak_width"]
    ic_rel_ht_filter = config_dict["ic_rel_height_filter"]
    ic_rel_ht_baseline = config_dict["ic_rel_height_filter_baseline"]
    ic_rel_ht_threshold = config_dict["ic_rel_height_threshold"]

    normalization_factors = pd.read_csv(args.normalization_factors_path)

    main(library_info_path=args.library_info_path,
         tensor_input_path=args.tensor_input_path,
         timepoints_dict=config_dict,
         normalization_factors=normalization_factors,
         isotope_clusters_out_path=args.isotope_clusters_out_path,
         factor_out_path=args.factor_data_out_path,
         factor_plot_output_path=args.factor_plot_out_path,
         ic_plot_output_path=args.ic_plot_out_path,
         gauss_params=args.gauss_params,
         filter_factors=filter_factors,
         factor_rt_r2_cutoff=factor_rt_r2_cutoff,
         factor_dt_r2_cutoff=factor_dt_r2_cutoff,
         ic_peak_prominence=ic_peak_prom,
         ic_peak_width=ic_peak_width,
         auto_ic_peak_width=auto_ic_peak_width,
         ic_rel_height_filter=ic_rel_ht_filter,
         ic_rel_height_filter_baseline=ic_rel_ht_baseline,
         ic_rel_height_threshold=ic_rel_ht_threshold)
