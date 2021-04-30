import os
import sys
import glob
import zlib
import math
import copy
import pickle
import pymzml
import molmass
import argparse
import importlib.util
import numpy as np
import pandas as pd
import seaborn as sns
import _pickle as cpickle
import matplotlib.pyplot as plt
from importlib import reload
from collections import defaultdict as ddict
from scipy.ndimage.filters import gaussian_filter

sys.path.append(os.getcwd() + "/workflow/scripts/")

from HDX_LIMIT.processing import TensorGenerator


#todo: later to use this function from hxtools(suggie version)
def calculate_theoretical_isotope_dist_from_sequence(sequence, n_isotopes=None):
    """Calculate theoretical isotope distribtuion from the given one-letter sequence of a library protein.

    Args:
        sequence (string): sequence in one letter code
        n_isotopes (int): number of isotopes to include. If none, includes all

    Return:
        isotope_dist (numpy ndarray): resulting theoretical isotope distribution

    """
    seq_formula = molmass.Formula(sequence)
    isotope_dist = np.array([x[1] for x in seq_formula.spectrum().values()])
    isotope_dist = isotope_dist / max(isotope_dist)
    if n_isotopes:
        if n_isotopes < len(isotope_dist):
            isotope_dist = isotope_dist[:n_isotopes]
        else:
            fill_arr = np.zeros(n_isotopes - len(isotope_dist))
            isotope_dist = np.append(isotope_dist, fill_arr)
    return isotope_dist


#todo: later to use this function from hxtools (suggie version)
def calculate_empirical_isotope_dist_from_integrated_mz(integrated_mz_array,
                                                        n_isotopes=None):
    """Calculate the isotope distribution from the integrated mz intensitities.
    
    Args: 
        integrated_mz_values (Numpy ndarray): array of integrated mz intensitites
        n_isotopes (int): number of isotopes to include. If none, includes all
    Returns: 
        isotope_dist (Numpy ndarray): isotope distribution with magnitude normalized to 1
    
    """
    isotope_dist = integrated_mz_array / max(integrated_mz_array)
    if n_isotopes:
        isotope_dist = isotope_dist[:n_isotopes]
    return isotope_dist


#todo: reminder! keep this function in this script rather than coming from hxtools
def calculate_isotope_dist_dot_product(sequence, undeut_integrated_mz_array):
    """Calculate dot product between theoretical isotope distribution from the sequence and experimental integrated mz array.
    
    Args:
        sequence (string): single-letter sequence of the library protein-of-interest
        undeut_integrated_mz_array (Numpy ndarray): observed integrated mz array from an undeuterated .mzML
    Returns:
        dot_product (float): result of dot product between theoretical and observed integrated-m/Z, from [0-1]

    """
    theo_isotope_dist = calculate_theoretical_isotope_dist_from_sequence(
        sequence=sequence)
    emp_isotope_dist = calculate_empirical_isotope_dist_from_integrated_mz(
        integrated_mz_array=undeut_integrated_mz_array)
    min_length = min([len(theo_isotope_dist), len(emp_isotope_dist)])
    dot_product = np.linalg.norm(
        np.dot(theo_isotope_dist[0:min_length], emp_isotope_dist[0:min_length])
    ) / np.linalg.norm(theo_isotope_dist) / np.linalg.norm(emp_isotope_dist)
    return dot_product


def gen_tensors_factorize(library_info_df,
                          undeut_tensor_path_list,
                          timepoint_index=0,
                          n_factors=15,
                          gauss_params=(3, 1)):
    """Instantiates TensorGenerator and factorizes.
    
    Args:
        library_info_df: library info data france
        undeut_tensor_path_list: undeuterated tensor file path list
        timepoint_index: time point index
        n_factors: number of factors for factorization
        gauss_params: gaussian paramters for factorization

    Return:
        undeut_ics_list (list): list of all IsotopeCluster objects from factorized tensors
        data_tensor_list (list): list of all DataTensor objects made from path_list

    """
    data_tensor_list = []
    undeut_ics_list = []

    for num, undeut_tensor_path in enumerate(undeut_tensor_path_list):

        # generate new data tensor
        new_data_tensor = TensorGenerator(filename=undeut_tensor_path,
                                          library_info=library_info_df,
                                          timepoint_index=timepoint_index)

        # factorize
        new_data_tensor.DataTensor.factorize(n_factors=n_factors,
                                             gauss_params=gauss_params)

        for factor in new_data_tensor.DataTensor.factors:
            for isotope_cluster in factor.isotope_clusters:
                undeut_ics_list.append(
                    isotope_cluster
                )

        data_tensor_list.append(new_data_tensor)

    return undeut_ics_list, data_tensor_list


def calc_dot_prod_for_isotope_clusters(sequence, undeut_isotope_clusters):
    """ 
    Calculate normalized dot product [0-1] of undeuterated IsotopeCluster.baseline_subtracted_int_mz to sequence determined theoretical distribution
    
    Parameters:
    sequence (str): sequence of the protein-of-interest in single-letter format
    undeut_isotope_clusters (list): list of IsotopeCluster objects to be compared against reference 

    Returns:
    dot_product_list (list): list of dot product results, index matched to integrated_mz_list
    integrated_mz_list (list): list of integrated m/Z arrays, index matched to dot_product_list

    """
    dot_product_list = []
    integrated_mz_list = []

    for index, isotope_clusters in enumerate(undeut_isotope_clusters):
        integrated_mz_array = isotope_clusters.baseline_integrated_mz
        dot_product = calculate_isotope_dist_dot_product(
            sequence=sequence, undeut_integrated_mz_array=integrated_mz_array)
        dot_product_list.append(dot_product)
        integrated_mz_list.append(integrated_mz_array)

    return dot_product_list, integrated_mz_list


def main(library_info_path,
         undeut_tensor_path_list,
         output_path=None,
         return_flag=None,
         n_factors=15,
         gauss_params=(3, 1)):
    """Compares each undeuterated replicate of a charge state to its theoretical distribution as a measure of signal quality.

    Args:
        library_info_path (string): path/to/library_info.csv
        undeut_tensor_path_list (list of strings): list of paths/to/files.cpickle.zlib
        output_path (string): path/to/output.csv
        return_flag (bool): option to return output in python, for notebook context
        n_factors (int): high number of factors to start factorization with
        gauss_params (tuple of floats): gaussian smoothing parameters in tuple (rt-sigma, dt-sigma), default (3,1)
    
    Returns:
        iso_cluster_list (list): list of IsotopeCluster objects produced from factorized input tensors
        data_tensor_list (list): list of DataTensor objects produced from input tensor paths
        idotp_list (list): list of resulting idotps for charge states
        integrated_mz_list (list): list containing integrated mzs of IsotopeClusters

    """
    print(undeut_tensor_path_list)
    lib_idx = int(undeut_tensor_path_list[0].split("/")[-1].split("_")[0])
    library_info = pd.read_csv(library_info_path)
    my_seq = library_info.iloc[lib_idx]["sequence"]

    iso_clusters_list, data_tensor_list = gen_tensors_factorize(
        library_info_df=library_info,
        undeut_tensor_path_list=undeut_tensor_path_list,
        n_factors=n_factors,
        gauss_params=gauss_params)

    idotp_list, integrated_mz_list = calc_dot_prod_for_isotope_clusters(
        sequence=my_seq, undeut_isotope_clusters=iso_clusters_list)

    if output_path is not None:
        pd.DataFrame({"idotp": max(idotp_list)}, index=[0]).to_csv(output_path)

    if return_flag is not None:
        return {
            "iso_clusters_list": iso_clusters_list,
            "data_tensor_list": data_tensor_list,
            "idotp_list": idotp_list,
            "integrated_mz_list": integrated_mz_list
        }


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=
        "Checks observed undeuterated signal against theoretical isotopic distribution, returns dataframe with highest idotp of all undeut"
    )
    parser.add_argument("library_info_path", help="path/to/library_info.csv")
    parser.add_argument(
        "undeut_tensor_path_list",
        nargs="+",
        help=
        "list of paths to undeuterated tensor outputs from extract_tensors.py")
    parser.add_argument("output_path", help="path/to/file for main .csv output")
    parser.add_argument(
        "--n_factors",
        default=15,
        help=
        "high number of factors to use in non_negative_parafac decomposition, counts down until correlation constraint is reached - see DataTensor.factorize()"
    )
    parser.add_argument("--gauss_params",
                        type=tuple,
                        default=(3, 1),
                        help="parameters for smoothing rt and dt dimensions")
    args = parser.parse_args()
    #### example of user inputs rather than from snakemake ####
    # library_info_path = '/Users/smd4193/Documents/MS_data/library_info.csv'
    # ins = ['/Users/smd4193/Documents/MS_data/1_20200922_lib15_2_0sec_01.mzML.gz.cpickle.zlib']
    # idotp_output = '/Users/smd4193/Documents/MS_data/1_idotp_check.csv'

    # main operation
    main(library_info_path=args.library_info_path,
         undeut_tensor_path_list=args.undeut_tensor_path_list,
         output_path=args.output_path)
