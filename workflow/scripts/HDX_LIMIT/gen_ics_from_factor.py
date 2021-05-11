import numpy as np
from scipy.signal import find_peaks
import _pickle as cpickle
import zlib
import copy
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error


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


    for factor in factor_data_dict['factors']:

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

        isotope_clusters['isotope_clusters'].append(ic_dict)


    print('heho')

    if output_path != None:
        save_ic_dict(obj=isotope_clusters, output_path=output_path)

    return isotope_clusters



def gauss_func(x, y0, A, xc, w):
    rxc = ((x - xc) ** 2) / (2 * (w ** 2))
    y = y0 + A * (np.exp(-rxc))
    return y


def estimate_gauss_param(array, xdata):
    ymax = np.max(array)
    maxindex = np.nonzero(array == ymax)[0]
    peakmax_x = xdata[maxindex][0]
    norm_arr = array/max(array)
    bins_for_width = norm_arr[norm_arr > 0.5]
    width_bin = len(bins_for_width)
    # binsnum = array[array > 0.2]
    # widthbin = len(binsnum[0])
    return peakmax_x, width_bin, ymax


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

    popt, pcov = curve_fit(gauss_func, xdata, ydata, method='lm', p0=[0, max_y, max_x, bin_width], maxfev=100000)

    y_fit = gauss_func(xdata, *popt)

    fit_mse = mean_squared_error(ydata, y_fit)

    gauss_fit_dict['data_label'] = data_label
    gauss_fit_dict['x_data'] = xdata
    gauss_fit_dict['y_data'] = ydata
    gauss_fit_dict['y_baseline'] = popt[0]
    gauss_fit_dict['y_amp'] = popt[1]
    gauss_fit_dict['xc'] = popt[2]
    gauss_fit_dict['width'] = popt[3]
    gauss_fit_dict['y_fit'] = y_fit
    gauss_fit_dict['fit_mse'] = fit_mse

    print('gauss_fit_dict')
    print(gauss_fit_dict)

    return gauss_fit_dict



def rel_height_peak_bounds(centers, int_mz, baseline_threshold=0.2, bound=30):
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



def find_isotope_cluster(factor_data_dict, int_threshold = 0.2, output_path=None):
    """
    find isotope cluster from factor data
    :param factors: factors dictionary from factor data dictionary
    :return:
    """

    # use integrated mz array to find peaks

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



        # creates ic from all factor data
        # if len(peaks) > 1:
        #     ic_idxs.append(
        #                 (feature_dict["left_bases"][0],
        #                  feature_dict["right_bases"][-1])
        #             )  # Create default ic from first left base to last right base

        # rel height bounds
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

    if output_path != None:
        save_ic_dict(isotope_clusters, output_path)

    return isotope_clusters


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

    factor_dict_fpath = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/factor.factor'
    factor_data_dict = load_factor_data(factor_dict_fpath)
    find_isotope_clusters_old(factor_data_dict, output_path=factor_dict_fpath+'oldic.ic')
    find_isotope_cluster(factor_data_dict, int_threshold=0.2, output_path=factor_dict_fpath+'.ic')
