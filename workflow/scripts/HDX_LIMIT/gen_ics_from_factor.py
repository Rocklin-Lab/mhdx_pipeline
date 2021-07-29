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
import glob
import numpy as np
from scipy.signal import find_peaks
import _pickle as cpickle
import zlib
import copy
import yaml
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from scipy.stats import linregress
import pandas as pd

from datatypes import IsotopeCluster
from processing import generate_tensor_factors, create_factor_data_object


def limit_write(obj, outpath):
    """ Writes a Python object as a zlib-compressed pickle.

	Args:
		obj (any Python Object):
		outpath (string): path/to/file.cpickle.zlib to write out

	Returns:
		None

	"""
    with open(outpath, "wb") as file:
        file.write(zlib.compress(cpickle.dumps(obj)))


def find_isotope_clusters_old(factor_data_dict, output_path=None):

    def rel_height_peak_bounds(centers, int_mz, bound=20):
        out = []
        baseline = max(int_mz) * 0.15  # TODO: HARDCODE
        for center in centers:
            if int_mz[center] > baseline:
                i, j = center, center
                cutoff = int_mz[center] * (bound / 100)
                while center - i <= 10 and i - 1 != -1:
                    i -= 1
                    if int_mz[i] < cutoff:
                        break
                while j - center <= 10 and j + 1 != len(int_mz):
                    j += 1
                    if int_mz[j] < cutoff:
                        break
                out.append((i, j))
        return out

    isotope_clusters = dict()
    isotope_clusters['mz_labels'] = factor_data_dict['mz_labels']
    isotope_clusters['bins_per_isotope_peak'] = factor_data_dict['bins_per_isotope_peak']
    isotope_clusters['isotope_clusters'] = []

    ics_idcs_list = []

    factors = factor_data_dict['factors']


    for factor in factors:

        ic_dict = dict()
        ic_dict['factor_mz_data'] = factor['factor_mz']
        integrated_mz = factor['factor_integrated_mz']

        ic_dict['factor_integrated_mz'] = integrated_mz

        peaks, feature_dict = find_peaks(integrated_mz,
                                     prominence=0.01,
                                     width=0.5)

        ic_dict['ic_mz_data_'] = []

        if len(peaks) == 0:
            return
        else:
            ic_idxs = [
                (feature_dict["left_bases"][i], feature_dict["right_bases"][i])
                for i in range(len(peaks))
                if
                feature_dict["left_bases"][i] < feature_dict["right_bases"][i]
                if feature_dict["right_bases"][i] -
                feature_dict["left_bases"][i] > 4
            ]
            # ic_idxs = [(feature_dict['left_bases'][i], feature_dict['left_bases'][i+1]) if feature_dict['left_bases'][i] < feature_dict['left_bases'][i+1] else (feature_dict['left_bases'][i], feature_dict['left_bases'][i]+6) for i in range(len(out[0])-1)]
            if len(peaks) > 1:
                ic_idxs.append(
                    (feature_dict["left_bases"][0],
                     feature_dict["right_bases"][-1])
                )  # Create default ic from first left base to last right base
            height_filtered = rel_height_peak_bounds(
                peaks, factor['factor_integrated_mz'])
            [ic_idxs.append(tup) for tup in height_filtered]
            cluster_idx = 0
            for integrated_indices in ic_idxs:
                if integrated_indices != None:
                    # try:

                    ic_mz_data = gen_isotope_peaks(factor['factor_mz'], integrated_indices,
                                                   factor_data_dict['bins_per_isotope_peak'])
                    isotope_peak_array = np.reshape(ic_mz_data, (-1, factor_data_dict['bins_per_isotope_peak']))

                    baseline_integrated_mz = np.sum(isotope_peak_array, axis=1)

                    peak_error = np.average(
                        np.abs(np.argmax(isotope_peak_array, axis=1) - ((factor_data_dict['bins_per_isotope_peak'] - 1) / 2)) / (
                                (factor_data_dict['bins_per_isotope_peak'] - 1) / 2), weights=baseline_integrated_mz)
                    baseline_peak_error = peak_error

                    outer_rtdt = sum(sum(np.outer(factor['factor_rt'], factor['factor_dt'])))

                    auc = sum(ic_mz_data) * outer_rtdt

                    if (baseline_peak_error / auc < 0.2):

                        ic_idx_dict = dict()
                        ic_dict['ic_mz_data_'].append(ic_mz_data)

                        ics_idcs_list.append(ic_idxs)

        isotope_clusters['isotope_clusters'].append(ic_dict)


    print('heho')

    if output_path != None:
        save_ic_dict(obj=isotope_clusters, output_path=output_path)

    return isotope_clusters, ics_idcs_list



def gauss_func(x, y0, A, xc, w):
    rxc = ((x - xc) ** 2) / (2 * (w ** 2))
    y = y0 + A * (np.exp(-rxc))
    return y

def gauss_func_no_baseline(x, A, xc, w):
    y0 = 0
    rxc = ((x - xc) ** 2) / (2 * (w ** 2))
    y = y0 + A * (np.exp(-rxc))
    return y


def estimate_gauss_param(array, xdata):
    ymax = np.max(array)
    maxindex = np.nonzero(array == ymax)[0]
    peakmax_x = xdata[maxindex][0]
    norm_arr = array/max(array)
    bins_for_width = norm_arr[norm_arr > 0.8]
    width_bin = len(bins_for_width)
    # binsnum = array[array > 0.2]
    # widthbin = len(binsnum[0])
    return peakmax_x, width_bin, ymax


def adjrsquared(r2, param, num):
    y = 1 - (((1 - r2) * (num - 1)) / (num - param - 1))
    return y


def fit_gaussian(xdata, ydata, data_label='dt'):

    gauss_fit_dict = dict()

    print('xdata')
    print(xdata)

    print('ydata')
    print(ydata)

    max_x, bin_width, max_y = estimate_gauss_param(ydata, xdata)

    print('max_x')
    print(max_x)

    print('bin_width')
    print(bin_width)

    print('max_y')
    print(max_y)

    popt, pcov = curve_fit(gauss_func, xdata, ydata, p0=[0, max_y, max_x, bin_width], maxfev=1000000)
    # popt, pcov = curve_fit(gauss_func_no_baseline, xdata, ydata, p0=[max_y, max_x, bin_width], maxfev=1000000)

    y_fit = gauss_func(xdata, *popt)
    # y_fit = gauss_func_no_baseline(xdata, *popt)

    fit_mse = mean_squared_error(ydata, y_fit)

    slope, intercept, rvalue, pvalue, stderr = linregress(ydata, y_fit)
    adj_r2 = adjrsquared(r2=rvalue**2, param=4, num=len(ydata))

    gauss_fit_dict['data_label'] = data_label
    gauss_fit_dict['x_data'] = xdata
    gauss_fit_dict['y_data'] = ydata
    gauss_fit_dict['y_baseline'] = 0
    gauss_fit_dict['y_baseline'] = popt[0]
    gauss_fit_dict['y_amp'] = popt[1]
    gauss_fit_dict['xc'] = popt[2]
    gauss_fit_dict['width'] = popt[3]
    gauss_fit_dict['y_fit'] = y_fit
    gauss_fit_dict['fit_mse'] = fit_mse
    gauss_fit_dict['fit_lingress_slope'] = slope
    gauss_fit_dict['fit_lingress_intercept'] = intercept
    gauss_fit_dict['fit_lingress_pvalue'] = pvalue
    gauss_fit_dict['fit_lingress_stderr'] = stderr
    gauss_fit_dict['fit_linregress_r2'] = rvalue**2
    gauss_fit_dict['fit_lingress_adj_r2'] = adj_r2

    print('gauss_fit_dict')
    print(gauss_fit_dict)

    return gauss_fit_dict



def rel_height_peak_bounds(centers, int_mz, baseline_threshold=0.2, bound=20):
        out = []
        baseline = max(int_mz) * baseline_threshold  # TODO: HARDCODE
        print('baseline: ', baseline)
        for center in centers:
            print('center: ', center)
            if int_mz[center] > baseline:
                i, j = center, center
                print('i: ', i)
                print('j: ', j)
                cutoff = int_mz[center] * (bound / 100)
                print('cutoff: ', cutoff)
                while center - i <= 10 and i - 1 != -1:
                    i -= 1
                    if int_mz[i] < cutoff:
                        break
                while j - center <= 10 and j + 1 != len(int_mz):
                    j += 1
                    if int_mz[j] < cutoff:
                        break
                out.append((i, j))
        return out


def gen_isotope_peaks(factor_mz_data, ic_idxs, bins_per_isotope_peak):

    print('ic_idxs')
    print(ic_idxs)

    low_idx = bins_per_isotope_peak * ic_idxs[0]
    high_idx = bins_per_isotope_peak * (ic_idxs[1] + 1)

    print('low_idx')
    print(low_idx)

    print('high_idx')
    print(high_idx)

    ic_mz_data = copy.deepcopy(factor_mz_data)
    ic_mz_data[0:low_idx] = 0
    ic_mz_data[high_idx:] = 0

    return ic_mz_data


def filter_factor_with_dt_rt_gauss_mse(factor_list, rt_mse_cutoff=0.1, dt_mse_cutoff=0.1):
    """
    curate factor based on rt and dt gauss fit
    :param factor_list:
    :param dt_mse_cutoff:
    :param rt_mse_cutoff:
    :return: curated factor list
    """

    filtered_factor_list = []

    for factor in factor_list:
        factor_rt_ind = np.arange(len(factor['factor_rt']))
        rt_gauss_fit = fit_gaussian(factor_rt_ind, factor['factor_rt'], data_label='rt')

        factor_dt_ind = np.arange(len(factor['factor_dt']))
        dt_gauss_fit = fit_gaussian(factor_dt_ind, factor['factor_dt'], data_label='dt')

        if rt_gauss_fit['fit_mse'] <= rt_mse_cutoff and dt_gauss_fit['fit_mse'] <= dt_mse_cutoff:

            filtered_factor_list.append(factor)

    new_factor_list = filtered_factor_list

    return new_factor_list



def filter_factor_with_dt_rt_gauss_r2(factor_list, rt_r2_cutoff=0.91, dt_r2_cutoff=0.91): # todo: see what a good value is here to filter
    """
    curate factor based on rt and dt gauss fit
    :param factor_list:
    :param dt_mse_cutoff:
    :param rt_mse_cutoff:
    :return: curated factor list
    """

    filtered_factor_list = []

    for factor in factor_list:
        factor_rt_ind = np.arange(len(factor['factor_rt']))
        rt_gauss_fit = fit_gaussian(factor_rt_ind, factor['factor_rt'], data_label='rt')

        factor_dt_ind = np.arange(len(factor['factor_dt']))
        dt_gauss_fit = fit_gaussian(factor_dt_ind, factor['factor_dt'], data_label='dt')

        if rt_gauss_fit['fit_linregress_r2'] >= rt_r2_cutoff:
            if dt_gauss_fit['fit_linregress_r2'] >= dt_r2_cutoff:

                filtered_factor_list.append(factor)

    new_factor_list = filtered_factor_list

    return new_factor_list



def make_new_dirpath(dirpath):

    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)

    return dirpath



def gen_isotope_cluster(tensor_fpath, library_info_fpath, timepoints_yaml, gauss_params, n_factors, factor_output_directory=None,
                        ic_output_directory=None, hdx_lim_ic_output_directory=None, filter_factors=False,
                        newICgen=True, oldICgen=True):

    factor_output_directory = make_new_dirpath(factor_output_directory)
    factor_plot_output_directory = make_new_dirpath(factor_output_directory+'/factor_plots')

    tensor_fname = os.path.split(tensor_fpath)[1]

    factor_output_fpath = os.path.join(factor_output_directory, tensor_fname+'.factor')
    factor_plot_output_path = os.path.join(factor_plot_output_directory, tensor_fname+'.factor.pdf')

    ic_output_directory = make_new_dirpath(ic_output_directory)

    hdx_lim_ic_output_directory = make_new_dirpath(hdx_lim_ic_output_directory)

    old_ic_dir = make_new_dirpath(hdx_lim_ic_output_directory + '/old_ics')
    old_ic_path = os.path.join(old_ic_dir, tensor_fname)

    new_ic_ht_filtered_dir = make_new_dirpath(hdx_lim_ic_output_directory + '/new_ics_ht_filt')
    new_ic_ht_filtered_path = os.path.join(new_ic_ht_filtered_dir, tensor_fname)

    new_ic_no_ht_filtered_dir = make_new_dirpath(hdx_lim_ic_output_directory + '/new_ics_no_ht_filt')
    new_ic_no_ht_filtered_path = os.path.join(new_ic_no_ht_filtered_dir, tensor_fname)

    if old_ic_path not in glob.glob(old_ic_dir + '/*.zlib') or new_ic_ht_filtered_path not in glob.glob(new_ic_ht_filtered_dir+'/*.zlib'):


        open_timepoints = yaml.load(open(timepoints_yaml, 'rb'), Loader=yaml.Loader)

        library_info = pd.read_csv(library_info_fpath)

        for tp in open_timepoints["timepoints"]:
            for fn in open_timepoints[tp]:
                if fn in tensor_fpath:
                    my_tp = tp

        print(my_tp)

        data_tensor = generate_tensor_factors(tensor_fpath=tensor_fpath,
                                              library_info_df=library_info,
                                              timepoint_index=my_tp,
                                              gauss_params=gauss_params,
                                              n_factors=n_factors,
                                              factor_output_fpath=None,
                                              factor_plot_output_path=None,
                                              timepoint_label=None)

        factor_data_dictionary = create_factor_data_object(data_tensor=data_tensor,
                                                           gauss_params=gauss_params,
                                                           timepoint_label=None)

        factors = factor_data_dictionary['factors']

        if filter_factors:
            # factors = filter_factor_with_dt_rt_gauss_mse(factors, dt_mse_cutoff=0.1, rt_mse_cutoff=0.1)
            factors = filter_factor_with_dt_rt_gauss_r2(factors, dt_r2_cutoff=0.91, rt_r2_cutoff=0.91)

        factor_data_dictionary['factors'] = factors


        if len(factors) > 0:


            if newICgen:

                ics_ht_filt_list = []
                ics_no_ht_filt_list = []

                new_IC_dir = make_new_dirpath(ic_output_directory + '/new_ICs')
                new_IC_fpath = os.path.join(new_IC_dir, tensor_fname + '.ic')

                iso_cluster_dict, ics_ht_filtered_idcs, ics_no_ht_filtered_idcs = find_isotope_cluster(factor_data_dictionary,
                                                                                                       output_path=new_IC_fpath)

                print('heho')

                iso_clusters = iso_cluster_dict['isotope_clusters']

                print('len iso_clusters')
                print(len(iso_clusters))

                print('len factors')
                print(len(data_tensor.DataTensor.factors))

                # todo: need to iterate over factor list coming from filter factor function below!!!!


                for num, factor in enumerate(factors):

                    ind_ht_filt_ics = ics_ht_filtered_idcs[num]
                    ind_no_ht_filt_ics = ics_no_ht_filtered_idcs[num]

                    for num2, ind_ht_filt in enumerate(ind_ht_filt_ics):

                        newIC = IsotopeCluster(
                            charge_states=data_tensor.DataTensor.charge_states,
                            factor_mz_data=copy.deepcopy(factor['factor_mz']),
                            source_file=data_tensor.DataTensor.source_file,
                            tensor_idx=data_tensor.DataTensor.tensor_idx,
                            timepoint_idx=data_tensor.DataTensor.timepoint_idx,
                            n_factors=n_factors,
                            factor_idx=factor['factor_num'],
                            cluster_idx=num2,
                            low_idx=data_tensor.DataTensor.bins_per_isotope_peak * ind_ht_filt[0],
                            high_idx=data_tensor.DataTensor.bins_per_isotope_peak * (ind_ht_filt[1] + 1),
                            rts=factor['factor_rt'],
                            dts=factor['factor_dt'],
                            retention_labels=data_tensor.DataTensor.retention_labels,
                            drift_labels=data_tensor.DataTensor.drift_labels,
                            mz_labels=data_tensor.DataTensor.mz_labels,
                            bins_per_isotope_peak=data_tensor.DataTensor.bins_per_isotope_peak,
                            max_rtdt=max(factor['factor_rt']) * max(factor['factor_rt']),
                            outer_rtdt=sum(sum(np.outer(factor['factor_rt'], factor['factor_rt']))),
                            n_concatenated=data_tensor.DataTensor.n_concatenated,
                            concat_dt_idxs=None,
                        )

                        ics_ht_filt_list.append(newIC)


                    for num3, ind_not_ht_filt in enumerate(ind_no_ht_filt_ics):

                        newIC = IsotopeCluster(
                            charge_states=data_tensor.DataTensor.charge_states,
                            factor_mz_data=copy.deepcopy(factor['factor_mz']),
                            source_file=data_tensor.DataTensor.source_file,
                            tensor_idx=data_tensor.DataTensor.tensor_idx,
                            timepoint_idx=data_tensor.DataTensor.timepoint_idx,
                            n_factors=n_factors,
                            factor_idx=factor['factor_num'],
                            cluster_idx=num3,
                            low_idx=data_tensor.DataTensor.bins_per_isotope_peak * ind_not_ht_filt[0],
                            high_idx=data_tensor.DataTensor.bins_per_isotope_peak * (ind_not_ht_filt[1] + 1),
                            rts=factor['factor_rt'],
                            dts=factor['factor_dt'],
                            retention_labels=data_tensor.DataTensor.retention_labels,
                            drift_labels=data_tensor.DataTensor.drift_labels,
                            mz_labels=data_tensor.DataTensor.mz_labels,
                            bins_per_isotope_peak=data_tensor.DataTensor.bins_per_isotope_peak,
                            max_rtdt=max(factor['factor_rt']) * max(factor['factor_rt']),
                            outer_rtdt=sum(sum(np.outer(factor['factor_rt'], factor['factor_rt']))),
                            n_concatenated=data_tensor.DataTensor.n_concatenated,
                            concat_dt_idxs=None,
                        )

                        ics_no_ht_filt_list.append(newIC)

            if oldICgen:

                old_ics_list = []

                old_IC_dir = make_new_dirpath(ic_output_directory + '/old_ICs')
                old_IC_fpath = os.path.join(old_IC_dir, tensor_fname + '.ic')

                iso_cluster_dict, ics_idcs = find_isotope_clusters_old(factor_data_dictionary, old_IC_fpath)

                for num, factor in enumerate(factors):

                    ind_ics = ics_idcs[num]

                    for num2, ind_ic in enumerate(ind_ics):

                        newIC = IsotopeCluster(
                            charge_states=data_tensor.DataTensor.charge_states,
                            factor_mz_data=copy.deepcopy(factor['factor_mz']),
                            source_file=data_tensor.DataTensor.source_file,
                            tensor_idx=data_tensor.DataTensor.tensor_idx,
                            timepoint_idx=data_tensor.DataTensor.timepoint_idx,
                            n_factors=n_factors,
                            factor_idx=factor['factor_num'],
                            cluster_idx=num2,
                            low_idx=data_tensor.DataTensor.bins_per_isotope_peak * ind_ic[0],
                            high_idx=data_tensor.DataTensor.bins_per_isotope_peak * (ind_ic[1] + 1),
                            rts=factor['factor_rt'],
                            dts=factor['factor_dt'],
                            retention_labels=data_tensor.DataTensor.retention_labels,
                            drift_labels=data_tensor.DataTensor.drift_labels,
                            mz_labels=data_tensor.DataTensor.mz_labels,
                            bins_per_isotope_peak=data_tensor.DataTensor.bins_per_isotope_peak,
                            max_rtdt=max(factor['factor_rt']) * max(factor['factor_rt']),
                            outer_rtdt=sum(sum(np.outer(factor['factor_rt'], factor['factor_rt']))),
                            n_concatenated=data_tensor.DataTensor.n_concatenated,
                            concat_dt_idxs=None,
                        )

                        old_ics_list.append(newIC)



            # save old ics
            # old_ic_dir = make_new_dirpath(hdx_lim_ic_output_directory + '/old_ics')
            # old_ic_path = os.path.join(old_ic_dir, tensor_fname)
            limit_write(old_ics_list, old_ic_path)

            # save new ics
            # ics with ht filtered
            # new_ic_ht_filtered_dir = make_new_dirpath(hdx_lim_ic_output_directory + '/new_ics_ht_filt')
            # new_ic_ht_filtered_path = os.path.join(new_ic_ht_filtered_dir, tensor_fname)
            limit_write(ics_ht_filt_list, new_ic_ht_filtered_path)

            # ics with no ht filtered
            # new_ic_no_ht_filtered_dir = make_new_dirpath(hdx_lim_ic_output_directory + '/new_ics_no_ht_filt')
            # new_ic_no_ht_filtered_path = os.path.join(new_ic_no_ht_filtered_dir, tensor_fname)
            limit_write(ics_no_ht_filt_list, new_ic_no_ht_filtered_path)





def find_isotope_cluster(factor_data_dict, int_threshold = 0.2, output_path=None):
    """
    find isotope cluster from factor data
    :param factors: factors dictionary from factor data dictionary
    :return:
    """

    # use integrated mz array to find peaks

    ics_ht_filtered_idcs = []
    ics_no_ht_filtered_idcs = []


    factors = factor_data_dict['factors']

    isotope_clusters = dict()
    isotope_clusters['mz_labels'] = factor_data_dict['mz_labels']
    isotope_clusters['bins_per_isotope_peak'] = factor_data_dict['bins_per_isotope_peak']
    isotope_clusters['isotope_clusters'] = []


    for factor in factors:

        ic_dict = dict()

        ic_dict['factor_mz_data'] = factor['factor_mz']

        factor_rt_ind = np.arange(len(factor['factor_rt']))
        rt_gauss_fit = fit_gaussian(factor_rt_ind, factor['factor_rt'], data_label='rt')



        factor_dt_ind = np.arange(len(factor['factor_dt']))
        dt_gauss_fit = fit_gaussian(factor_dt_ind, factor['factor_dt'], data_label='dt')


        ic_dict['factor_rt_gauss_fit'] = rt_gauss_fit
        ic_dict['factor_dt_gauss_fit'] = dt_gauss_fit

        integrated_mz = factor['factor_integrated_mz']
        norm_integrated_mz = integrated_mz/max(integrated_mz)

        ic_dict['factor_integrated_mz'] = integrated_mz

        peaks, feature_dict = find_peaks(norm_integrated_mz, prominence=int_threshold)

        ic_idxs = [
                        (feature_dict["left_bases"][i], feature_dict["right_bases"][i])
                        for i in range(len(peaks))
                        if
                        feature_dict["left_bases"][i] < feature_dict["right_bases"][i]
                        if feature_dict["right_bases"][i] -
                        feature_dict["left_bases"][i] > 4
                    ]


        height_filtered = rel_height_peak_bounds(
                    peaks, norm_integrated_mz)

        ic_dict['ic_mz_data_ht_filtered'] = []
        for num, ht_filt_idx in enumerate(height_filtered):
            ic_mz_data = gen_isotope_peaks(factor_mz_data=factor['factor_mz'],
                                           ic_idxs=ht_filt_idx,
                                           bins_per_isotope_peak=factor_data_dict['bins_per_isotope_peak'])
            ic_dict['ic_mz_data_ht_filtered'].append(ic_mz_data)


        ic_dict['ic_mz_data_no_ht_filtered'] = []

        for ind, ic_idx_no_ht in enumerate(ic_idxs):
            ic_mz_data = gen_isotope_peaks(factor_mz_data=factor['factor_mz'],
                                           ic_idxs=ic_idx_no_ht,
                                           bins_per_isotope_peak=factor_data_dict['bins_per_isotope_peak'])
            ic_dict['ic_mz_data_no_ht_filtered'].append(ic_mz_data)

        isotope_clusters['isotope_clusters'].append(ic_dict)

        ics_ht_filtered_idcs.append(height_filtered)
        ics_no_ht_filtered_idcs.append(ic_idxs)

    if output_path != None:
        save_ic_dict(isotope_clusters, output_path)

    return isotope_clusters, ics_ht_filtered_idcs, ics_no_ht_filtered_idcs


def save_ic_dict(obj, output_path):

    with open(output_path, "wb") as file:
        file.write(zlib.compress(cpickle.dumps(obj)))



def load_factor_data(factor_data_filepath):
    """
    plot factor data from factor data file .factor
    :param factor_data_filepath: .factor filepath
    :return: None. saves the figure
    """


    factor_data = cpickle.loads(zlib.decompress(open(factor_data_filepath, 'rb').read()))

    return factor_data


if __name__ == '__main__':

    # NO HARDCODED PATHS! THESE VALUES SHOULD BE HANDLED BY ARGPARSE!
    # TODO: Set up argparse for CLI inputs.

    factor_dict_fpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/factor.factor'
    tensor_fpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/tensor_files/old_tensor_files/_eehee_rd4_0871_tensor_files/113_20210323_lib15_ph6_0sec_01.mzML.gz.cpickle.zlib'
    library_info_fpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/library_info/old_library_info/library_info.csv'
    config_fpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/config.yaml'
    # factor_data_dict = load_factor_data(factor_dict_fpath)
    # find_isotope_clusters_old(factor_data_dict, output_path=factor_dict_fpath+'oldic.ic')
    # find_isotope_cluster(factor_data_dict, int_threshold=0.2, output_path=factor_dict_fpath+'_filter_factor.ic')

    tensor_files_dpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/tensor_files/_eehee_rd4_0871_tensor_files'
    tensor_files_list = glob.glob(tensor_files_dpath + '/*.gz.cpickle.zlib')

    factor_output_dpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/tensor_files/_eehee_rd4_0871_tensor_files/gauss_filter/factors'
    ic_output_dpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/tensor_files/_eehee_rd4_0871_tensor_files/gauss_filter/ics_suggie'
    hdx_lim_ic_dpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/tensor_files/_eehee_rd4_0871_tensor_files/gauss_filter/ics_hdxlim'

    for ind, tensor_fpath in enumerate(tensor_files_list):

        gen_isotope_cluster(tensor_fpath=tensor_fpath,
                            library_info_fpath=library_info_fpath,
                            timepoints_yaml=config_fpath,
                            gauss_params=(3, 1),
                            n_factors=15,
                            factor_output_directory=factor_output_dpath,
                            ic_output_directory=ic_output_dpath,
                            hdx_lim_ic_output_directory=hdx_lim_ic_dpath,
                            filter_factors=True,
                            newICgen=True,
                            oldICgen=True)
