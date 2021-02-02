# do all the hx rate fitting to the data here

import importlib.util

auxiliary_spec = importlib.util.spec_from_file_location(
    "auxiliary.py", "./auxiliary/rate_fitting/auxiliary.py"
)
hx_rate_fit_functions_spec = importlib.util.spec_from_file_location(
    "molmass.py", "./auxiliary/rate_fitting/hx_rate_fit_functions.py"
)
auxiliary = importlib.util.module_from_spec(auxiliary_spec)
hx_rate_fit_functions = importlib.util.module_from_spec(hx_rate_fit_functions_spec)
auxiliary_spec.loader.exec_module(auxiliary)
hx_rate_fit_functions_spec.loader.exec_module(hx_rate_fit_functions)

import os
import numpy as np
from scipy.special import expit
from sklearn.metrics import mean_squared_error
from scipy.optimize import fmin_powell, basinhopping
from auxiliary import (
    calc_rmse,
    load_pickle_file,
    calculate_theoretical_isotope_dist_from_sequence,
)
from hx_rate_fit_functions import (
    isotope_dist_from_PoiBin_calc_with_hx_rates,
    calculate_intrinsic_exchange_rates_suggie,
)
import Bio.PDB
import matplotlib.pyplot as plt


def calc_back_exchange_without_free_energies(
    theo_isotope_dist, exp_isotope_dist, intrinsic_rates, init_guess=95
):
    """
    calculate the back exchange. model back exchange according to the error distribution
    :param theo_isotope_dist:
    :param exp_isotope_dist:
    :param intrinsic_rates:
    :return:
    """
    opt = fmin_powell(
        lambda x: calc_rmse(
            isotope_dist_from_PoiBin_calc_with_hx_rates(
                isotope_dist=theo_isotope_dist,
                timepoints=1e9,
                rates=intrinsic_rates,
                num_bins=len(exp_isotope_dist),
                backexchange=expit(x),
            ),
            exp_isotope_dist,
        ),
        x0=init_guess,
        disp=False,
    )
    back_exchange = expit(opt)

    return back_exchange


def poibin_isotope_dist_concat(
    timepoints, rates, num_bins, exp_isotope_dist_concat, backexchange_arr
):
    """
    calculate the poibin isotope dist
    :param timepoints: timepoints
    :param rates: intrinsic rates
    :param num_bins: number of bins to include for isotope dist
    :param exp_isotope_dist_concat: experimental isotope dist concatenated
    :param backexchange: backexchange
    :return: poibin isotope dist
    """
    out_arr = np.zeros((len(timepoints), num_bins))
    for ind, timepoint in enumerate(timepoints):
        out_arr[ind, :] = isotope_dist_from_PoiBin_calc_with_hx_rates(
            exp_isotope_dist_concat,
            timepoint,
            rates,
            num_bins=num_bins,
            backexchange=backexchange_arr[ind],
        )
    out_arr = np.ravel(out_arr)
    return out_arr


def hx_rate_fit_rmse(
    timepoints,
    rates,
    thr_isotope_dist_list,
    exp_isotope_dist_concat,
    num_bins,
    backexchange_arr,
):
    """

    :param timepoints: timepoints
    :param rates: rates
    :param thr_isotope_dist_list: theoretical isotope dist list
    :param num_bins: number of bins
    :param backexchange: backexchange
    :return: rmse between exp isotope dist concat and concat model from PoiBin Isotope dist
    """
    concat_model = poibin_isotope_dist_concat(
        timepoints, rates, num_bins, thr_isotope_dist_list, backexchange_arr
    )
    mean_sq_error = mean_squared_error(
        exp_isotope_dist_concat[exp_isotope_dist_concat > 0],
        concat_model[exp_isotope_dist_concat > 0],
    )
    return mean_sq_error


def fit_hx_rates_optimize(
    timepoints,
    thr_isotope_dist_list,
    exp_isotope_dist_list,
    num_rates,
    backexchange,
    num_bins=None,
    n_iter=200,
    temp=0.00003,
    step_size=0.02,
):
    """
    calculate the hx rates using the time points and the list of experimental isotope distribution
    :param timepoints: timepoints
    :param exp_isotope_dist_list: list of experimental isotope distributions
    :param backexchange: backexchange rate
    :return: fitted hx rates
    """

    backexchange_arr = np.array([backexchange for x in timepoints])
    if len(backexchange_arr) != len(timepoints):
        backexchange_arr = np.reshape(backexchange_arr, (len(timepoints),))

    exp_isotope_dist_concat = np.concatenate(exp_isotope_dist_list)

    init_rates = np.linspace(-7, 0, num_rates)

    opt = basinhopping(
        lambda rates: hx_rate_fit_rmse(
            timepoints=timepoints,
            rates=rates,
            thr_isotope_dist_list=thr_isotope_dist_list,
            exp_isotope_dist_concat=exp_isotope_dist_concat,
            num_bins=num_bins,
            backexchange_arr=backexchange_arr,
        ),
        x0=init_rates,
        niter=n_iter,
        T=temp,
        stepsize=step_size,
        minimizer_kwargs={"options": {"maxiter": 1}},
    )
    return opt


def plot_rates(intrinsic_rate, fitted_rate, output_path, pdb_name, extension=".png"):
    """
    plot intrinsic and fitted rate
    :param intrinsic_rate: sorted intrinsic rate list/array
    :param fitted_rate: sorted fitted rate list/array
    :param output_path: directory to save figure
    :param pdb_name: pdb name
    :param extension: file extension (png, pdf)
    :return: None. Saves the output file
    """
    plt.plot(intrinsic_rate, color="red", label="intrinsic_rate")
    plt.plot(fitted_rate, color="blue", label="fitted_rate")
    plt.ylabel("rate")
    plt.xlabel("res")
    plt.legend(loc="best", fontsize="small")
    plt.savefig(os.path.join(output_path, pdb_name + "_hx_rate" + extension))
    plt.close()


def write_hx_rate_output(
    intrinsic_rates,
    fitted_rates,
    output_path,
    pdb_name,
    pD,
    temp,
    backexchange,
    temp_bh,
    step_size_bh,
    n_iter_bh,
):
    """
    write the hx rate output to a csv file
    :param intrinsic_rate: sorted intrinsic rate list/array
    :param fitted_rate: sorted fitted rate list/array
    :param output_path: dirpath to save output
    :param pdb_name: pdb name
    :param pD: pD
    :param temp: temperature
    :param temp_bh: temp for basin hopping
    :param step_size_bh: step size for basin opping
    :param n_iter_bh: number of iterations for basin hopping
    :return: None. Saves the output to .csv file
    """
    param_main_header = "PARAMS\n"
    param_header = "pdb_name,pD,temp,backexchange,temp_bh,n_iter_bh,step_size_bh\n"
    param_data = "{},{},{},{},{},{},{}\n".format(
        pdb_name,
        str(pD),
        str(temp),
        str(backexchange),
        str(temp_bh),
        str(n_iter_bh),
        str(step_size_bh),
    )

    # data section
    data_main_header = "\nOUTPUT\n"
    data_header = "num,intrinsic_rates,fitted_rates\n"
    data_string = ""

    for ind, (intrin_rate, fitted_rate) in enumerate(
        zip(intrinsic_rates, fitted_rates)
    ):
        line = "{},{},{}\n".format(str(ind + 1), intrin_rate, fitted_rate)
        data_string += line

    output_string = (
        param_main_header
        + param_header
        + param_data
        + data_main_header
        + data_header
        + data_string
    )

    with open(os.path.join(output_path, pdb_name + "_hx_rate.csv"), "w") as outfile:
        outfile.write(output_string)
        outfile.close()


def fit_hx_rates(
    pdb_path,
    timepoint_list,
    ms_data_list,
    output_path,
    num_bins=50,
    pD=6.15,
    temp=298,
    nterm=None,
    cterm=None,
    temp_bh=0.0003,
    step_size_bh=0.02,
    n_iter_bh=200,
    write_output=False,
    plot_rate=False,
):
    """
    fit hx rates to individual residues, saves outputs, and returns the fitted rate
    :param pdb_path: path to pdb structure
    :param timepoint_list: timepoint list
    :param ms_data_list: ms data list
    :param output_path: dirpath to output files
    :param num_bins: number of bins to include in the isotope distribution
    :param pD: pD
    :param temp: temperature in Kelvin
    :param nterm: any addition to N-terminus of the structure
    :param cterm: any addition to C-terminus of the structure
    :param temp_bh: temperature set for basin hopping hx rate fitting optimization
    :param step_size_bh: step size set for basin hopping hx rate fitting optimization
    :param n_iter_bh: number of iterations for basin hopping hx rate fitting optimization
    :param write_output: write the hx rate to a .csv file
    :param plot_rate: boolean. If True, will output plot of intrinsic vs fitted rate
    :return: sorted fitted rate
    """
    # load the pdb structure using BioPython and generate a string of one letter residue sequence
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path, pdb_path)
    residues = [res for res in structure[0].get_residues()]
    sequence = "".join(
        [Bio.PDB.Polypeptide.three_to_one(aa.get_resname()) for aa in residues]
    )

    # add any additional residues at nterminal or cterminal ends. These variables are important later on while calculating
    # structural parameters such as burial depth, secondary structure elements, etc.

    if nterm:
        sequence = nterm + sequence
    if cterm:
        sequence = sequence + cterm

    # calculate the theoretical isotope distribution from the given sequence
    theoretical_isotope_dist = calculate_theoretical_isotope_dist_from_sequence(
        sequence
    )

    # calculate intrinsic rates based on the sequence, temperature, and pD
    intrinsic_rates = calculate_intrinsic_exchange_rates_suggie(
        sequence, Temperature=temp, pH=pD
    )
    intrinsic_rates[1] = 0.0  # can't properly calculate

    # normalize the ms data and make a list of isotope distribution
    isotope_dist_list = []
    for ind, ms_data in enumerate(ms_data_list):
        isotope_dist_list.append(ms_data / max(ms_data))

    # calculate back exchange rate using the last (or the fully deuterated**) isotopde distribution
    back_exchange = calc_back_exchange_without_free_energies(
        theo_isotope_dist=theoretical_isotope_dist,
        exp_isotope_dist=isotope_dist_list[-1],
        intrinsic_rates=intrinsic_rates,
        init_guess=2,
    )

    # calculate how many rates to use for fitting (non-zero ones)
    len_rates = len([x for x in intrinsic_rates if x != 0])

    # fitting rate using basin hopping optimization
    fit_rates = fit_hx_rates_optimize(
        timepoints=timepoint_list,
        thr_isotope_dist_list=theoretical_isotope_dist,
        exp_isotope_dist_list=isotope_dist_list,
        num_rates=len_rates,
        backexchange=back_exchange,
        num_bins=num_bins,
        n_iter=n_iter_bh,
        temp=temp_bh,
        step_size=step_size_bh,
    )

    # sort the fitted rates in an ascending manner
    sorted_fitted_hx_rates = sorted(fit_rates.x)

    # writes the output to a csv file if True
    if write_output:
        sorted_intrinsic_rates = sorted(np.log(intrinsic_rates[intrinsic_rates > 0]))
        pdb_path_components = os.path.split(pdb_path)
        pdb_name = pdb_path_components[1]

        write_hx_rate_output(
            sorted_intrinsic_rates,
            sorted_fitted_hx_rates,
            output_path,
            pdb_name,
            pD,
            temp,
            back_exchange,
            temp_bh,
            step_size_bh,
            n_iter_bh,
        )

    # plots the figure if True
    if plot_rate:
        pdb_path_components = os.path.split(pdb_path)
        pdb_name = pdb_path_components[1]
        sorted_intrinsic_rates = sorted(np.log(intrinsic_rates[intrinsic_rates > 0]))
        plot_rates(
            sorted_intrinsic_rates, sorted_fitted_hx_rates, output_path, pdb_name
        )

    return sorted_fitted_hx_rates


if __name__ == "__main__":

    dirpath = "/Users/smd4193/Documents/hx_ratefit_gabe"
    pickle_path = os.path.join(dirpath, "fp3_n_1717_HHH_rd4_0997.pdb_z5.0.pickle")
    pdb_path = os.path.join(dirpath, "HHH_rd4_0997.pdb")

    timepoints = [
        8,
        9,
        13,
        20,
        30,
        60,
        100,
        180,
        310,
        480,
        13 * 60,
        19 * 60,
        34 * 60,
        57 * 60,
        82 * 60,
        158 * 60,
        268 * 60,
        ((7 * 60) + 49) * 60,
        ((13 * 60) + 24) * 60,
        ((18 * 60) + 24) * 60,
        ((26 * 60) + 50) * 60,
    ]  # in seconds

    pickle_obj = load_pickle_file(pickle_path)
    pickle_ms_data = pickle_obj[1]

    ms_data_list = []
    for ind in range(len(pickle_ms_data)):
        ms_data_list.append(pickle_ms_data[ind]["comp_filtered"])

    fitted_hx_rate_sorted = fit_hx_rates(
        pdb_path=pdb_path,
        timepoint_list=timepoints,
        ms_data_list=ms_data_list[1:],
        output_path=dirpath,
        num_bins=50,
        pD=6.15,
        temp=298,
        nterm="HM",
        cterm="GSSGSSGNS",
        temp_bh=0.00003,
        step_size_bh=0.02,
        n_iter_bh=20,
        write_output=True,
        plot_rate=False,
    )

    ###############################################
    ## used for testing
    # dirpath = "/Users/smd4193/Documents/hx_ratefit_gabe"
    # pickle_path = os.path.join(dirpath, "fp3_n_1717_HHH_rd4_0997.pdb_z5.0.pickle")
    # pdb_path = os.path.join(dirpath, "HHH_rd4_0997.pdb")
    # pdb_patah_20053 = os.path.join(dirpath, "HHH_rd4_0518")
    # Nterm = 'HM'
    # Cterm = 'GSSGSSGNS'
    # structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path, pdb_path)
    # residues = [res for res in structure[0].get_residues()]
    # sequence = ''.join([Bio.PDB.Polypeptide.three_to_one(aa.get_resname()) for aa in residues])
    # sequence = Nterm + sequence + Cterm
    #
    # pD = 6.15
    # temp = 295
    # intrinsic_rates = calculate_intrinsic_exchange_rates_suggie(sequence, Temperature=temp, pH=pD)
    # intrinsic_rates[1] = 0.0
    #
    # timepoints = [8, 9, 13, 20, 30, 60,
    #               100, 180, 310, 480,
    #               13 * 60, 19 * 60, 34 * 60, 57 * 60,
    #               82 * 60, 158 * 60, 268 * 60, ((7 * 60) + 49) * 60,
    #               ((13 * 60) + 24) * 60, ((18 * 60) + 24) * 60, ((26 * 60) + 50) * 60]  # in seconds
    #
    # pickle_obj = load_pickle_file(pickle_path)
    #
    # datafile = dict()
    #
    # datafile['output'] = pickle_obj[1]
    # isotope_dist_list = []
    # for ind in range(len(datafile['output'])):
    #     datafile['output'][ind]['major_species_integrated_intensities'] = datafile['output'][ind]['comp_filtered']
    #     isotope_dist_list.append(datafile['output'][ind]['major_species_integrated_intensities'] / max(
    #         datafile['output'][ind]['major_species_integrated_intensities']))
    #
    # exp_isotope_dist_final_timepoint = datafile['output'][-1]['major_species_integrated_intensities']/max(datafile['output'][-1]['major_species_integrated_intensities'])
    #
    # theo_isotope = calculate_theoretical_isotope_dist_from_sequence(sequence)
    # # back_exch = 0.88
    # back_exch = calc_back_exchange_without_free_energies(theo_isotope,
    #                                                      exp_isotope_dist_final_timepoint,
    #                                                      intrinsic_rates=intrinsic_rates,
    #                                                      init_guess=2)
    #
    # len_rates = len([x for x in intrinsic_rates if x!=0])
    #
    # fit_rates = fit_hx_rates_optimize(timepoints, theo_isotope, isotope_dist_list[1:],
    #                          num_rates=len_rates,
    #                          backexchange=back_exch,
    #                          n_iter=5,
    #                          num_bins=50)
    #
    # print('backexch ', back_exch)
    # print(fit_rates.x)
    # plt.plot(sorted(fit_rates.x), color='blue')
    # plt.plot(sorted(np.log(intrinsic_rates[intrinsic_rates > 0])), color='red')
    # plt.show()
    # plt.close()
    # print('gheho')
    #
    # print('heho')
