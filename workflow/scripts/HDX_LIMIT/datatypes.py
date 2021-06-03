import os
import time
import sys
import math
import copy
import psutil
import peakutils
import numpy as np
from nn_fac import ntf
from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.measurements import center_of_mass
import scipy as sp
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from scipy.stats import norm


class DataTensor:

    def __init__(self, source_file, tensor_idx, timepoint_idx, name,
                 total_mass_window, n_concatenated, charge_states, integrated_mz_limits, bins_per_isotope_peak, normalization_factor, **kwargs):

        ###Set Common Attributes###

        # Positional args
        self.source_file = source_file
        self.tensor_idx = tensor_idx
        self.timepoint_idx = timepoint_idx
        self.name = name
        self.total_mass_window = total_mass_window
        self.n_concatenated = n_concatenated
        self.charge_states = charge_states
        self.integrated_mz_limits = integrated_mz_limits
        self.bins_per_isotope_peak = bins_per_isotope_peak
        self.normalization_factor = normalization_factor

        # Keyword Args
        if kwargs is not None:
            kws = list(kwargs.keys())
            if "rts" in kws:
                self.rts = np.array(kwargs["rts"])
            if "dts" in kws:
                self.dts = np.array(kwargs["dts"])
            if "seq_out" in kws:
                self.seq_out = np.array(kwargs["seq_out"])
            if "int_seq_out" in kws and kwargs["int_seq_out"] is not None:
                self.int_seq_out = np.array(kwargs["int_seq_out"])
                self.int_seq_out_float = self.int_seq_out.astype("float64")
                self.int_grid_out = np.reshape(
                    self.int_seq_out_float, (len(self.rts), len(self.dts), 50))
                self.int_gauss_grids = self.gauss(self.int_grid_out)
            if "concat_dt_idxs" in kws:
                self.concat_dt_idxs = kwargs["concat_dt_idxs"]
            if "concatenated_grid" in kws:
                self.concatenated_grid = kwargs["concatenated_grid"]


        ###Compute Instance Values###

        # Handle normal case: new DataTensor from output of isolate_tensors.py
        if self.n_concatenated == 1:

            (
                self.retention_labels,
                self.drift_labels,
                self.mz_labels,
                self.bins_per_isotope_peak,
                self.full_grid_out,
            ) = self.sparse_to_full_tensor_reprofile((self.rts, self.dts, self.seq_out),
                                                     self.integrated_mz_limits,
                                                     self.bins_per_isotope_peak )
            self.full_gauss_grids = self.gauss(self.full_grid_out)

        # Handle concatenated tensor case, check for required inputs
        else:
            if not all("dts" in kwargs and "rts" in kwargs and
                       "lows" in kwargs and "highs" in kwargs and
                       "concatenated_grid" in kwargs and
                       "abs_mz_low" in kwargs and "concat_dt_idxs" in kwargs):

                print("Concatenated Tensor Missing Required Values")
                sys.exit()

    def sparse_to_full_tensor_reprofile(self, data, integrated_mz_limits, bins_per_isotope_peak = 7, ms_resolution=25000):
        retention_labels, drift_labels, sparse_data = data
        
        FWHM = np.average(integrated_mz_limits) / ms_resolution
        gaussian_scale = FWHM / 2.355 #a gaussian with scale = standard deviation = 1 has FWHM 2.355
        
        mz_bin_centers = np.ravel([np.linspace(lowlim, highlim, bins_per_isotope_peak) for lowlim, highlim in integrated_mz_limits])
        tensor3_out = np.zeros((len(retention_labels), len(drift_labels), len(mz_bin_centers)))


        scan = 0
        for i in range(len(retention_labels)):
            for j in range(len(drift_labels)):
                n_peaks = len(sparse_data[scan])
                gaussians = norm(loc=sparse_data[scan][:,0],scale=gaussian_scale)
                resize_gridpoints = np.resize(mz_bin_centers, (n_peaks, len(mz_bin_centers) )).T
                eval_gaussians = gaussians.pdf(resize_gridpoints) * sparse_data[scan][:,1] * gaussian_scale

                tensor3_out[i][j] = np.sum(eval_gaussians,axis=1) 
                scan += 1
                
        return retention_labels, drift_labels, mz_bin_centers, bins_per_isotope_peak, tensor3_out


    # Takes tensor input and gaussian filter parameters, outputs filtered data
    def gauss(self, grid, rt_sig=3, dt_sig=1):

        gauss_grid = np.zeros(np.shape(grid))
        for i in range(np.shape(grid)[2]):
            gauss_grid[:, :, i] = gaussian_filter(grid[:, :, i],
                                                  (rt_sig, dt_sig))
        return gauss_grid

    def factorize(self, n_factors=4, new_mz_len=None, gauss_params=None): 
        # Test factorization starting at n_factors = 15 and counting down, keep factorization that has no factors with correlation greater than 0.2 in any dimension.

        def corr_check(factors):
            # Checks scipy non_negatve_parafac output factors for inter-factor (off-diagonal) correlations > cutoff, returns True if all values are < cutoff

            a = np.minimum(
                np.minimum(np.corrcoef(factors[0].T),
                           np.corrcoef(factors[1].T)),
                np.corrcoef(factors[2].T),
            )

            return np.max(a[np.where(~np.eye(a.shape[0], dtype=bool))])

        def pmem(id_str):
            process = psutil.Process(os.getpid())
            print(id_str + " Process Memory (GB): " +
                  str(process.memory_info().rss / 1024 / 1024 / 1024))

        t = time.time()
        pmem("0 Start")
        # print('Filtering... T+'+str(t-t0))
        # handle concatenation and intetrpolfilter option
        if self.n_concatenated != 1:
            #code handing n_concatenated != 1 needs  to be re-written from scratch
            grid, lows, highs, concat_dt_idxs = (
                self.concatenated_grid,
                self.concat_dt_idxs,
            )
        else:
            concat_dt_idxs = None
            if gauss_params != None:
                grid = self.gauss(self.full_grid_out, gauss_params[0],
                                  gauss_params[1])
            else:
                grid = self.full_grid_out
        
        grid = self.full_gauss_grids
        
        pmem("1 Pre-Factorization")
        n_itr = 2
        
        last_corr_check = 1.0
        n_factors += 1
        while n_factors > 2 and last_corr_check > 0.17:
            n_factors -= 1
            pmem(str(n_itr) + " " + str(n_factors) + " Factors " + " Start")
            t1 = time.time()
            # print('Starting '+str(nf)+' Factors... T+'+str(t1-t))
            nnf1 = ntf.ntf(grid, n_factors)
            pmem(str(n_itr) + " " + str(n_factors) + " Factors " + " End")
            n_itr += 1
            t2 = time.time()
            # print('Factorization Duration: '+str(t2-t1))

            if n_factors > 1:
                last_corr_check = corr_check(nnf1)
                
        pmem(str(n_itr) + " Post-Factorization")
        n_itr += 1
        # Create Factor objects
        factors = []
        t = time.time()
        # print('Saving Factor Objects... T+'+str(t-t0))
        for i in range(n_factors):
            pmem(str(n_itr) + " Start Factor " + str(i))
            n_itr += 1
            factors.append(
                Factor(
                    source_file=self.source_file,
                    tensor_idx=self.tensor_idx,
                    timepoint_idx=self.timepoint_idx,
                    name=self.name,
                    charge_states=self.charge_states,
                    rts=nnf1[0].T[i],
                    dts=nnf1[1].T[i],
                    mz_data=nnf1[2].T[i],
                    retention_labels=self.retention_labels,
                    drift_labels=self.drift_labels,
                    mz_labels=self.mz_labels,
                    factor_idx=i,
                    n_factors=n_factors,
                    bins_per_isotope_peak = self.bins_per_isotope_peak,
                    n_concatenated=self.n_concatenated,
                    concat_dt_idxs=concat_dt_idxs,
                    normalization_factor=self.normalization_factor
                ))
            pmem(str(n_itr) + " End Factor " + str(i))
            n_itr += 1
        pmem(str(n_itr) + " Factor Initialization End")
        n_itr += 1
        self.factors = factors
        pmem(str(n_itr) + " Script End")
        # t = time.time()
        # print('Done: T+'+str(t-t0))



def cal_area_under_curve_from_normal_distribution(low_bound, upper_bound, center, width):
    """
    calculate area under the curve given the lower and upper bound
    :param low_bound: low bound
    :param upper_bound: upper bound
    :param center: center of distribution
    :param width: width of distribution
    :return: area under curve
    """

    lb_cdf = norm.cdf(low_bound, loc=center, scale=width)
    ub_cdf = norm.cdf(upper_bound, loc=center, scale=width)
    auc = ub_cdf - lb_cdf
    return auc


def estimate_gauss_param(ydata, xdata):
    ymax = np.max(ydata)
    maxindex = np.nonzero(ydata == ymax)[0]
    peakmax_x = xdata[maxindex][0]
    norm_arr = ydata/max(ydata)
    bins_for_width = norm_arr[norm_arr > 0.8]
    width_bin = len(bins_for_width)
    init_guess = [0, ymax, peakmax_x, width_bin]
    return init_guess


def gauss_func(x, y0, A, xc, w):
    rxc = ((x - xc) ** 2) / (2 * (w ** 2))
    y = y0 + A * (np.exp(-rxc))
    return y


def adjrsquared(r2, param, num):
    y = 1 - (((1 - r2) * (num - 1)) / (num - param - 1))
    return y


def fit_gaussian(xdata, ydata, data_label='dt'):

    init_guess = estimate_gauss_param(ydata, xdata)

    gauss_fit_dict = dict()
    gauss_fit_dict['data_label'] = data_label

    try:
        popt, pcov = curve_fit(gauss_func, xdata, ydata, p0=init_guess, maxfev=1000000)
        y_fit = gauss_func(xdata, *popt)
        fit_rmse = mean_squared_error(ydata/max(ydata), y_fit/max(y_fit), squared=False)
        slope, intercept, rvalue, pvalue, stderr = linregress(ydata, y_fit)
        adj_r2 = adjrsquared(r2=rvalue**2, param=4, num=len(ydata))
        gauss_fit_dict['gauss_fit_success'] = True
        gauss_fit_dict['y_baseline'] = popt[0]
        gauss_fit_dict['y_amp'] = popt[1]
        gauss_fit_dict['xc'] = popt[2]
        gauss_fit_dict['width'] = popt[3]
        gauss_fit_dict['y_fit'] = y_fit
        gauss_fit_dict['fit_rmse'] = fit_rmse
        gauss_fit_dict['fit_lingress_slope'] = slope
        gauss_fit_dict['fit_lingress_intercept'] = intercept
        gauss_fit_dict['fit_lingress_pvalue'] = pvalue
        gauss_fit_dict['fit_lingress_stderr'] = stderr
        gauss_fit_dict['fit_linregress_r2'] = rvalue ** 2
        gauss_fit_dict['fit_lingress_adj_r2'] = adj_r2
    except:
        gauss_fit_dict['gauss_fit_success'] = False
        gauss_fit_dict['xc'] = None
        gauss_fit_dict['width'] = None
        gauss_fit_dict['fit_rmse'] = 100
        gauss_fit_dict['fit_linregress_r2'] = 0.0

    return gauss_fit_dict


###
### Class - Factor:
### Holds data from a single component of a non_negative_parafac factorization
### constructed as: (nnf1[0].T[i], nnf1[1].T[i], nnf1[2].T[i], i, n, self.lows, self.highs, self.n_concatenated)
###


class Factor:

    def __init__(
        self,
        source_file,
        tensor_idx,
        timepoint_idx,
        name,
        charge_states,
        rts,
        dts,
        mz_data,
        retention_labels,
        drift_labels,
        mz_labels,
        factor_idx,
        n_factors,
        bins_per_isotope_peak,
        n_concatenated,
        concat_dt_idxs,
        normalization_factor
    ):

        ###Set Attributes###

        self.source_file = source_file
        self.tensor_idx = tensor_idx
        self.timepoint_idx = timepoint_idx
        self.name = name
        self.charge_states = charge_states
        self.rts = rts
        self.dts = dts
        self.mz_data = mz_data
        self.retention_labels = retention_labels
        self.drift_labels = drift_labels
        self.mz_labels = mz_labels
        self.auc = sum(mz_data)
        self.factor_idx = factor_idx
        self.n_factors = n_factors
        self.bins_per_isotope_peak = bins_per_isotope_peak
        self.n_concatenated = n_concatenated
        self.concat_dt_idxs = concat_dt_idxs
        self.normalization_factor = normalization_factor

        ###Compute Instance Values###

        # integrate within expected peak bounds and create boolean mask of expected peak bounds called grate
        self.integrated_mz_data = np.sum(np.reshape(mz_data, (-1, self.bins_per_isotope_peak)), axis=1)

        self.max_rtdt = max(self.rts) * max(self.dts)
        self.outer_rtdt = sum(sum(np.outer(self.rts, self.dts)))

        # This can be a shared function
        #self.integrated_mz_baseline = peakutils.baseline(
        #    np.asarray(self.integrated_mz_data),
        #    6)  # 6 degree curve seems to work well
        #       
        #self.baseline_subtracted_integrated_mz = (self.integrated_mz_data -
        #                                          self.integrated_mz_baseline)
        self.baseline_subtracted_integrated_mz = self.integrated_mz_data

        # fit factor rts and dts to gaussian
        self.rt_gauss_fit = fit_gaussian(np.arange(len(self.rts)), self.rts, data_label='rt')
        self.dt_gauss_fit = fit_gaussian(np.arange(len(self.dts)), self.dts, data_label='dt')


        # calculate rt and dt auc
        if self.rt_gauss_fit['gauss_fit_success']:
            self.rt_auc = cal_area_under_curve_from_normal_distribution(low_bound=0,
                                                                               upper_bound=len(self.rts) - 1,
                                                                               center=self.rt_gauss_fit['xc'],
                                                                               width=self.rt_gauss_fit['width'])
        else:
            self.rt_auc = None

        if self.dt_gauss_fit['gauss_fit_success']:
            self.dt_auc = cal_area_under_curve_from_normal_distribution(low_bound=0,
                                                                               upper_bound=len(self.dts) - 1,
                                                                               center=self.dt_gauss_fit['xc'],
                                                                               width=self.dt_gauss_fit['width'])
        else:
            self.rt_auc = None

        ## old protocol.
        # Writes to self.isotope_clusters
        # self.find_isotope_clusters()  # heuristic height value, should be high-level param TODO - Will require passage through DataTensor class

        ## now generates self.isotope_clusters upon calling the function self.find_isotope_clusters


    def rel_height_peak_bounds(self, centers, norm_integrated_mz, baseline_threshold=0.15, rel_ht_threshold=0.2):
        out = []
        for center in centers:
            if norm_integrated_mz[center] > baseline_threshold:
                i, j = center, center
                cutoff = norm_integrated_mz[center] * rel_ht_threshold
                while center - i <= 10 and i - 1 != -1:
                    i -= 1
                    if norm_integrated_mz[i] < cutoff:
                        break
                while j - center <= 10 and j + 1 != len(norm_integrated_mz):
                    j += 1
                    if norm_integrated_mz[j] < cutoff:
                        break
                out.append((i, j))
        return out



    # Uses find_window function to identify portions of the integrated mz dimension that look 'isotope-cluster-like', saves as Factor attribute
    def find_isotope_clusters(self, prominence=0.15, width_val=3, rel_height_filter=True, baseline_threshold=0.15, rel_height_threshold=0.10):
        """
        find isotope clusters by finding peaks in integrated mz distribution and choosing indices to include in mz data for ic clusters
        :param prominence: prominence value for choosing peaks.
        :param width_val: minimum width for choosing peaks
        :param rel_height_filter: within peaks, whether or not
        :param baseline_threshold: baseline for rel height filtering
        :param rel_height_threshold: rel height threhold for rel height filtering
        :return:
        """

        self.isotope_clusters = []

        norm_integrated_mz = self.baseline_subtracted_integrated_mz/max(self.baseline_subtracted_integrated_mz)

        peaks, feature_dict = find_peaks(norm_integrated_mz,
                                         prominence=prominence,
                                         width=width_val)

        if len(peaks) == 0:
            ic_idxs = [(0, len(self.baseline_subtracted_integrated_mz)-1)]
            int_mz_width = [2]
            # return
        else:
            int_mz_width = [
                feature_dict['widths'][i]
                for i in range(len(peaks))
                if
                feature_dict["left_bases"][i] < feature_dict["right_bases"][i]
                if feature_dict["right_bases"][i] -
                feature_dict["left_bases"][i] > 4
            ]
            ic_idxs = [
                (feature_dict["left_bases"][i], feature_dict["right_bases"][i])
                for i in range(len(peaks))
                if
                feature_dict["left_bases"][i] < feature_dict["right_bases"][i]
                if feature_dict["right_bases"][i] -
                feature_dict["left_bases"][i] > 4
            ]

            if rel_height_filter:
                height_filtered = self.rel_height_peak_bounds(centers=peaks,
                                                              norm_integrated_mz=norm_integrated_mz,
                                                              baseline_threshold=baseline_threshold,
                                                              rel_ht_threshold=rel_height_threshold)

                ic_idxs = height_filtered

            # ic_idxs = [(feature_dict['left_bases'][i], feature_dict['left_bases'][i+1]) if feature_dict['left_bases'][i] < feature_dict['left_bases'][i+1] else (feature_dict['left_bases'][i], feature_dict['left_bases'][i]+6) for i in range(len(out[0])-1)]
            ## previous protocol to create ic indexes to include the entire factor mz data. Not including this in ic generation.
            # if len(peaks) > 1:
            #     ic_idxs.append(
            #         (feature_dict["left_bases"][0],
            #          feature_dict["right_bases"][-1])
            #     )  # Create default ic from first left base to last right base
            # height_filtered = rel_height_peak_bounds(
            #     peaks, self.baseline_subtracted_integrated_mz)
            # [ic_idxs.append(tup) for tup in height_filtered]

        cluster_idx = 0
        for integrated_indices, integrated_mz_width in zip(ic_idxs, int_mz_width):
            if integrated_indices != None:
                #try:

                newIC = IsotopeCluster(
                    integrated_mz_peak_width=integrated_mz_width,
                    charge_states=self.charge_states,
                    factor_mz_data=copy.deepcopy(self.mz_data),
                    source_file=self.source_file,
                    tensor_idx=self.tensor_idx,
                    timepoint_idx=self.timepoint_idx,
                    n_factors=self.n_factors,
                    factor_idx=self.factor_idx,
                    cluster_idx=cluster_idx,
                    low_idx = self.bins_per_isotope_peak * integrated_indices[0],
                    high_idx = self.bins_per_isotope_peak * (integrated_indices[1] + 1),
                    rts=self.rts,
                    dts=self.dts,
                    retention_labels=self.retention_labels,
                    drift_labels=self.drift_labels,
                    mz_labels=self.mz_labels,
                    bins_per_isotope_peak=self.bins_per_isotope_peak,
                    max_rtdt=self.max_rtdt,
                    outer_rtdt=self.outer_rtdt,
                    n_concatenated=self.n_concatenated,
                    concat_dt_idxs=self.concat_dt_idxs,
                    normalization_factor=self.normalization_factor
                )
                if (newIC.baseline_peak_error / newIC.baseline_auc <
                        0.2):  # TODO: HARDCODE
                    self.isotope_clusters.append(newIC)
                    cluster_idx += 1
                #except:
                #print("ic index out of bounds: " + str(integrated_indices))
        return

    # heuristically identifies 'things that look like acceptable isotope clusters' in integrated mz dimension, roughly gaussian allowing some inflection points from noise
    def find_window(self, array, peak_idx, width):
        rflag = True
        lflag = True
        if peak_idx == 0:
            win_low = 0
            lflag = False
        if peak_idx == len(array) - 1:
            win_high = len(array) - 1
            rflag = False

        idx = peak_idx + 1
        if idx < len(array) - 1:  # Check if idx is last idx
            if (array[idx] < array[peak_idx] /
                    5):  # Peak is likely not an IC if peak > 5 x neighbors
                win_high = idx
                rflag = False

        while rflag:
            # make sure looking ahead won't throw error
            if idx + 1 < len(array):
                # if idx+1 goes down, and its height is greater than 20% of the max peak
                if array[idx + 1] < array[idx] and array[
                        idx + 1] > array[peak_idx] / 5:
                    idx += 1
                # if above check fails, test conditions separately
                else:
                    if array[idx + 1] < array[peak_idx] / 5:
                        win_high = idx
                        rflag = False
                    else:
                        # first check if upward point is more than 5x the height of the base peak
                        if array[idx + 1] < array[peak_idx] * 5:
                            # look one point past upward point, if its below the last point and the next point continues down, continue
                            if idx + 2 < len(array):
                                if array[idx + 2] < array[idx + 1]:
                                    if idx + 3 < len(array):
                                        if array[idx + 3] < array[idx + 2]:
                                            idx += 3
                                        else:  # point 3 ahead goes up, do not keep sawtooth tail
                                            win_high = idx
                                            rflag = False
                                    else:  # point 2 past idx is end of array, set as high limit
                                        win_high = idx + 2
                                        rflag = False
                                else:  # points continue increasing two ahead of idx, stop at idx
                                    win_high = idx
                                    rflag = False
                            else:  # upward point is end of array, do not keep
                                win_high = idx
                                rflag = False
                        else:  # upward point is major spike / base peak is minor, end search
                            win_high = idx
                            rflag = False
            else:  # idx is downward and end of array
                win_high = idx
                rflag = False

        idx = peak_idx - 1
        if idx >= 0:
            if array[idx] < array[peak_idx] / 5:
                win_low = idx
                lflag = False

        while lflag:
            if idx - 1 >= 0:  # make sure looking ahead won't throw error
                if (
                        array[idx - 1] < array[idx] and
                        array[idx - 1] > array[peak_idx] / 5
                ):  # if idx-1 goes down, and its height is greater than 20% of the max peak
                    idx -= 1
                # if above check fails, test conditions separately
                else:
                    if array[idx - 1] < array[peak_idx] / 5:
                        win_low = idx
                        lflag = False
                    else:
                        # first check if upward point is more than 5x the height of the base peak
                        if array[idx - 1] < array[peak_idx] * 5:
                            # look one point past upward point, if its below the last point and the next point continues down, continue
                            if idx - 2 >= 0:
                                if array[idx - 2] < array[idx]:
                                    if idx - 3 >= 0:
                                        if array[idx - 3] < array[idx - 2]:
                                            idx -= 3
                                        else:  # point 3 ahead goes up, do not keep sawtooth tail
                                            win_low = idx
                                            lflag = False
                                    else:  # point 2 past idx is end of array, set as high limit
                                        win_low = idx - 2
                                        lflag = False
                                else:  # points continue increasing two ahead of idx, stop at idx
                                    win_low = idx
                                    lflag = False
                            else:  # upward point is start of array, do not keep
                                win_low = idx
                                lflag = False
                        else:  # upward point is major spike / base peak is minor, end search
                            win_low = idx
                            lflag = False
            else:  # idx is start of array
                win_low = idx
                lflag = False

        if win_high - win_low < width:
            return None
        else:
            return [win_low, win_high]


###
### Class - IsotopeCluster:
### Contains mz data of mz region identified to have isotope-cluster-like characteristics, stores full data of parent factor
###


class IsotopeCluster:

    def __init__(
        self,
        integrated_mz_peak_width,
        charge_states,
        factor_mz_data,
        source_file,
        tensor_idx,
        timepoint_idx,
        n_factors,
        factor_idx,
        cluster_idx,
        low_idx,
        high_idx,
        rts,
        dts,
        retention_labels,
        drift_labels,
        mz_labels,
        bins_per_isotope_peak,
        max_rtdt,
        outer_rtdt,
        n_concatenated,
        concat_dt_idxs,
        normalization_factor
    ):

        ###Set Attributes###
        self.integrated_mz_peak_width = integrated_mz_peak_width
        self.charge_states = charge_states
        self.factor_mz_data = factor_mz_data
        self.source_file = source_file
        self.tensor_idx = tensor_idx
        self.timepoint_idx = timepoint_idx
        self.n_factors = n_factors
        self.factor_idx = factor_idx
        self.cluster_idx = cluster_idx
        self.low_idx = low_idx
        self.high_idx = high_idx
        self.rts = rts
        self.dts = dts
        self.retention_labels = retention_labels
        self.drift_labels = drift_labels
        self.mz_labels = mz_labels
        self.bins_per_isotope_peak = bins_per_isotope_peak
        self.max_rtdt = max_rtdt
        self.outer_rtdt = outer_rtdt
        self.n_concatenated = n_concatenated
        self.concat_dt_idxs = concat_dt_idxs
        self.normalization_factor = normalization_factor

        ###Calculate Scoring Requirements###

        # prune factor_mz to get window around cluster that is consistent between charge-states
        self.cluster_mz_data = copy.deepcopy(self.factor_mz_data)
        self.cluster_mz_data[0:self.low_idx] = 0
        self.cluster_mz_data[self.high_idx:] = 0

        # integrate area of IC
        self.auc = sum(self.cluster_mz_data) * self.outer_rtdt / self.normalization_factor

        # identify peaks and find error from expected peak positions using raw mz
        
        #CONTINUE HERE TO DEFINE PEAK ERROR USING bins_per_isotope_peak

        isotope_peak_array = np.reshape(self.cluster_mz_data, (-1, self.bins_per_isotope_peak))

        self.baseline = 0
        self.baseline_subtracted_mz = self.cluster_mz_data 
        self.baseline_auc = self.auc
        self.log_baseline_auc = np.log(self.baseline_auc)
        self.baseline_max_peak_height = max(self.baseline_subtracted_mz)
        self.baseline_integrated_mz = np.sum(isotope_peak_array, axis=1)

        self.peak_error = np.average( np.abs( np.argmax(isotope_peak_array,axis=1) - ((self.bins_per_isotope_peak - 1)/2) ) / ((self.bins_per_isotope_peak - 1)/2), weights=self.baseline_integrated_mz)
        self.baseline_peak_error = self.peak_error


        # Cache int_mz and rt scoring values
        # Compute baseline_integrated_mz_com, baseline_integrated_mz_std, baseline_integrated_mz_FWHM and baseline_integrated_mz_rmse from Gaussian Fit
        # Define functions for Gaussian fit BEGIN
        def gaussian_function(x, H, A, x0, sigma):
            return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

        def gauss_fit(x, y):
            mean = sum(x * y) / sum(y)
            sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
            nonzeros = [index for index, value in enumerate(list(y)) if value!=0]
            popt, pcov = curve_fit(gaussian_function, x, y, p0=[0, max(y), mean, sigma], bounds=([0, 0, nonzeros[0], 0], [np.inf, np.inf, nonzeros[-1], np.inf]))
            return popt

        def params_from_gaussian_fit(self):
            try:
                xdata = [i for i in range(len(self.baseline_integrated_mz))]
                ydata = self.baseline_integrated_mz
                H, A, x0, sigma = gauss_fit(xdata, ydata)
                y_gaussian_fit = gaussian_function(xdata, *gauss_fit(xdata, ydata))
                rmse = mean_squared_error(ydata/max(ydata), y_gaussian_fit/max(y_gaussian_fit), squared=False)
                com = x0
                std = sigma
                FWHM = 2*np.sqrt(2*np.log(2))*std
                return rmse, com, std, FWHM
            except:
                return 100, 100, 100, 100
            # Define functions for Gaussian fit END       
        
        
        self.baseline_integrated_mz_norm = self.baseline_integrated_mz / np.linalg.norm(
            self.baseline_integrated_mz)
        #self.baseline_integrated_mz_com = center_of_mass(
        #    self.baseline_integrated_mz)[
        #        0]  # COM in IC integrated bin dimension
        #self.baseline_integrated_mz_std = (np.average(
        #    (np.arange(len(self.baseline_integrated_mz)) -
        #     self.baseline_integrated_mz_com)**2,
        #    weights=self.baseline_integrated_mz,
        #)**0.5)
        self.baseline_integrated_mz_rmse, self.baseline_integrated_mz_com, self.baseline_integrated_mz_std, self.baseline_integrated_mz_FWHM = params_from_gaussian_fit(self)

        self.rt_norm = self.rts / np.linalg.norm(self.rts)
        self.rt_com = center_of_mass(self.rts)[0]

        # Cache DT values
        # If DT is concatenated, return list of coms and norms of single rts relative to bin numbers, a single_dt distribution starts at 0. If only one charge state, return list of len=1
        if self.concat_dt_idxs is not None:
            single_dts = []
            # generate list of single dts
            single_dts.append(self.dts[:self.concat_dt_idxs[0]])
            for i in range(len(self.charge_states) - 1):
                single_dts.append(
                    self.dts[self.concat_dt_idxs[i]:self.concat_dt_idxs[i + 1]])

            self.single_dts = single_dts
            self.dt_coms = [center_of_mass(dt)[0] for dt in single_dts]
            self.dt_norms = [dt / np.linalg.norm(dt) for dt in single_dts]
        else:
            self.dt_coms = [center_of_mass(self.dts)[0]]
            self.dt_norms = [self.dts / np.linalg.norm(self.dts)]

        if self.n_concatenated == 1:
            self.abs_mz_com = np.average(self.mz_labels, weights=self.cluster_mz_data)
        else:
            self.abs_mz_com = "Concatenated, N/A, see IC.baseline_integrated_mz_com"

        # format useful values to be read by pandas
        self.info_tuple = (
            self.
            source_file,  # Filename of data used to create parent DataTensor
            self.tensor_idx,  # Library master list row of parent-DataTensor
            self.n_factors,  # Number of factors in parent decomposition
            self.
            factor_idx,  # Index of IC parent-factor in DataTensor.factors[]
            self.cluster_idx,  # Index of IC in parent-factor.isotope_clusters[]
            self.charge_states,  # List of charge states in IC
            self.
            n_concatenated,  # number of source tensors IC parent-DataTensor was made from
            self.low_idx,  # Low bin index corresponding to Factor-level bins
            self.high_idx,  # High bin index corresponding to Factor-level bins
            self.baseline_auc,  # Baseline-subtracted AUC (BAUC)
            self.baseline_auc,  # Baseline-subtracted grate area sum (BGS) #gabe 210507: duplicating baseline_auc for now bc we got rid of grate sum
            self.
            baseline_peak_error,  # Baseline-subtracted version of peak-error (BPE)
            self.
            baseline_integrated_mz_com,  # Center of mass in added mass units
            self.
            abs_mz_com,  # IsotopicCluster center-of-mass in absolute m/z dimension
            self.rts,  # Array of RT values
            self.
            dts,  # Array of DT values, if a tensor is concatenated,this is taken from the last tensor in the list, can be seen in tensor_idx
            np.arange(0, len(self.baseline_integrated_mz), 1),
            self.
            baseline_integrated_mz,  # Array of baseline-subtracted integrated mz intensity values
        )

        # instantiate to make available for setting in PathOptimizer
        self.bokeh_tuple = None  # bokeh plot info tuple
        self.single_sub_scores = None  # differences between winning score and score if this IC were substituted, list of signed values
        self.undeut_ground_dot_products = None
