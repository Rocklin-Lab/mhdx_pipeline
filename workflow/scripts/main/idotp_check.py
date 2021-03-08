import os
import sys
import glob
import zlib
import math
import copy
import pickle
import pymzml
import importlib.util
import numpy as np
import pandas as pd
import seaborn as sns
import _pickle as cpickle
import matplotlib.pyplot as plt
from importlib import reload
from collections import defaultdict as ddict
from scipy.ndimage.filters import gaussian_filter

sys.path.append(os.getcwd()+"/workflow/scripts/auxiliary/")

#the following sys path is for local debug
# sys.path.append('../auxiliary/')

import LC_IM_MS_TensorAnalysis as hx
import molmass


#todo: later to use this function from hxtools(suggie version)
def calculate_theoretical_isotope_dist_from_sequence(sequence, num_isotopes=None):
    """
    calculate theoretical isotope distribtuion from a given one letter sequence of protein chain
    :param sequence: sequence in one letter code
    :param num_isotopes: number of isotopes to include. If none, includes all
    :return: isotope distribution
    """
    seq_formula = molmass.Formula(sequence)
    isotope_dist = np.array([x[1] for x in seq_formula.spectrum().values()])
    isotope_dist = isotope_dist/max(isotope_dist)
    if num_isotopes:
        if num_isotopes < len(isotope_dist):
            isotope_dist = isotope_dist[:num_isotopes]
        else:
            fill_arr = np.zeros(num_isotopes - len(isotope_dist))
            isotope_dist = np.append(isotope_dist, fill_arr)
    return isotope_dist

#todo: later to use this function from hxtools (suggie version)
def calculate_empirical_isotope_dist_from_integrated_mz(integrated_mz_array, num_isotopes=None):
    """
    calculate the isotope distribution from the integrated mz intensitities
    :param integrated_mz_values: array of integrated mz intensitites
    :param num_isotopes: number of isotopes to include. If none, includes all
    :return: isotope distribution
    """
    isotope_dist = integrated_mz_array/max(integrated_mz_array)
    if num_isotopes:
        isotope_dist = isotope_dist[:num_isotopes]
    return isotope_dist


#todo: reminder! keep this function in this script rather than coming from hxtools
def calculate_isotope_dist_dot_product(sequence, undeut_integrated_mz_array):
    """
    calculate dot product between theoretical isotope distribution from the sequence and experimental integrated mz array
    :param sequence: sequence of the protein
    :param undeut_integrated_mz_array: integrated mz array
    :return: dot product
    """
    theo_isotope_dist = calculate_theoretical_isotope_dist_from_sequence(sequence=sequence)
    emp_isotope_dist = calculate_empirical_isotope_dist_from_integrated_mz(integrated_mz_array=undeut_integrated_mz_array)
    min_length = min([len(theo_isotope_dist), len(emp_isotope_dist)])
    dot_product = np.linalg.norm(np.dot(theo_isotope_dist[0:min_length], emp_isotope_dist[0:min_length])) / np.linalg.norm(theo_isotope_dist) / np.linalg.norm(emp_isotope_dist)
    return dot_product


def gen_tensors_factorize(library_info_df, undeut_tensor_fpath_list, timepoint_index=0, num_factors=13, factor_gauss_param=(3, 1)):
    """
    generate data tensor and factorizes
    :param library_info_df: library info data france
    :param undeut_tensor_fpath_list: undeuterated tensor file path list
    :param timepoint_index: time point index
    :param num_factors: number of factors for factorization
    :param factor_gauss_param: gaussian paramters for factorization
    :return: list of isotope_cluster and data_tensor
    """

    gauss_undeut_ics_list = []
    data_tensor_list = []

    for num, undeut_tensor_fpath in enumerate(undeut_tensor_fpath_list):

        # generate new data tensor
        new_data_tensor = hx.TensorGenerator(filename=undeut_tensor_fpath,
                                          library_info=library_info_df,
                                          timepoint_index=timepoint_index)

        # factorize
        new_data_tensor.DataTensor.factorize(n_factors=num_factors,
                                          gauss_params=factor_gauss_param)

        for factor in new_data_tensor.DataTensor.factors:
            for isotope_cluster in factor.isotope_clusters:
                gauss_undeut_ics_list.append(isotope_cluster) #todo: check if the outer bracket creates a different list shape
        
        data_tensor_list.append(new_data_tensor)

    return gauss_undeut_ics_list, data_tensor_list


def calc_dot_prod_for_isotope_clusters(sequence, gauss_undeut_isotope_clusters):
    """
    calculate dot products for all the factorized isotope clusters
    :param sequence: sequence of the protein
    :param gauss_undeut_isotope_clusters: isotope cluster from factorization
    :return: dot product list between comparision of theoretical isotope distribution and isotope clusters from
    factorization + list of integrated mz list
    """

    dot_product_list = []
    integrated_mz_list = []

    for index, isotope_clusters in enumerate(gauss_undeut_isotope_clusters):
        integrated_mz_array = isotope_clusters.baseline_integrated_mz
        dot_product = calculate_isotope_dist_dot_product(sequence=sequence,
                                                      undeut_integrated_mz_array=integrated_mz_array)
        dot_product_list.append(dot_product)
        integrated_mz_list.append(integrated_mz_array)

    return dot_product_list, integrated_mz_list


def gen_idotp_check_dataframe(library_info_fpath, undeut_tensor_fpath_list, output_fpath, timepoint_index=0, num_factors=13,
                              factor_gauss_param=(3,1)):
    """

    :param library_info_fpath: library info file path
    :param undeut_tensor_fpath_list: undeut tensor filepath list
    :param output_fpath: idotp check output file path
    :param timepoint_index: timepoint index
    :param num_factors: number of factors for factorization
    :param factor_gauss_param: gauss param for factorization
    :return: iso_cluster_list, data_tensor_list, idotp_list, integrated_mz_list (for debugging / checking)
    """

    lib_idx = int(undeut_tensor_fpath_list[0].split('/')[-1].split('_')[0])
    library_info = pd.read_csv(library_info_fpath)
    my_seq = library_info.iloc[lib_idx]['sequence']

    iso_clusters_list, data_tensor_list = gen_tensors_factorize(library_info_df=library_info,
                                                                undeut_tensor_fpath_list=undeut_tensor_fpath_list,
                                                                timepoint_index=timepoint_index,
                                                                num_factors=num_factors,
                                                                factor_gauss_param=factor_gauss_param)

    idotp_list, integrated_mz_list = calc_dot_prod_for_isotope_clusters(sequence=my_seq,
                                                                        gauss_undeut_isotope_clusters=iso_clusters_list)

    pd.DataFrame({'idotp': max(idotp_list)}, index=[0]).to_csv(output_fpath)

    return iso_clusters_list, data_tensor_list, idotp_list, integrated_mz_list


if __name__ == '__main__':

    # input and output paths determined here. You can replace the snakemake inputs to put your own filepaths
    library_info_fpath = snakemake.input.pop(0)
    ins = snakemake.input
    idotp_output = snakemake.output[0]

    #### example of user inputs rather than from snakemake ####
    # library_info_fpath = '/Users/smd4193/Documents/MS_data/library_info.csv'
    # ins = ['/Users/smd4193/Documents/MS_data/1_20200922_lib15_2_0sec_01.mzML.gz.cpickle.zlib']
    # idotp_output = '/Users/smd4193/Documents/MS_data/1_idotp_check.csv'

    # main operation
    idotp_out = gen_idotp_check_dataframe(library_info_fpath=library_info_fpath,
                              undeut_tensor_fpath_list=ins,
                              output_fpath=idotp_output,
                              timepoint_index=0,
                              num_factors=13,
                              factor_gauss_param=(3,1))
