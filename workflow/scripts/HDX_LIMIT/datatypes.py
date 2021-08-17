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
from scipy.stats import linregress


class DataTensor:
    """A container for LC-IMS-MS data that includes smoothing and factorization methods.

    Attributes:
        source_file (str): Path of DataTensor's parent resources/tensors/.cpickle.zlib file.
        tensor_idx (int): Index of the DataTensor in a concatenated DataTensor (concatenation of tensors is deprecated, this value always 0).
        timepoint_idx (int): Index of the DataTensors HDX timepoint in config["timepoints"].
        name (str): Name of DataTensor's rt-group.
        total_mass_window (int): Magnitude of DataTensor's m/Z dimension in bins.
        n_concatenated (int): Number of DataTensors combined to make the instance DataTensor. Deprecated, value always 1.
        charge_states (list with one int): Net positive charge on protein represented in DataTensor, format left over from concatenation.
        integrated_mz_limits (numpy 2D array): Values for low and high m/Z limits of integration around expected peak centers.
        bins_per_isotope_peak (int): Number of integrated m/Z bins for a signal to be considered as an IsotopeCluster.
        normalization_factor (float): Divisor for the integrated m/Z intensity of any IsotopeCluster from a parent DataTensor instance. 
        rts (numpy array): Intensity summed over each rt bin.
        dts (numpy array): Intensity summed over each dt bin.
        seq_out (numpy array): Intensity summed over each m/Z bin in a flat array.
        int_seq_out (numpy array): Intensities integrated over all dimensions in a flat array. (TODO: Review, not sure about this)
        int_seq_out_float (numpy array of float64):  Intensities integrated over all dimensions in a flat array, cast as float64.
        int_grid_out (numpy array): int_seq_out reshaped to a 3D array by rt, dt, and m/Z dimension magnitudes.
        int_gauss_grids (numpy array): int_grid_out after applying gaussian smoothing.
        concat_dt_idxs (list of ints): Deprecated - Indices marking the boundaries between different DataTensors' concatenated DT dimensions.
        concatenated_grid (numpy array): Deprecated - int_grid_out from multiple DataTensors concatenated along the dt axis.
        retention_labels (list of floats): Mapping of DataTensor's RT bins to corresponding absolute retention time in minutes.
        drift_labels (list of floats): Mapping of DataTensor's DT bins to corresponding absolute drift time in miliseconds.
        mz_labels (list of floats): Mapping of DataTensor's m/Z bins to corresponding m/Z.
        full_grid_out (numpy array): Reprofiled tensor with intensities integrated within bounds defined in integrated_mz_limits.
        full_gauss_grids (numpy array): full_grid_out after the application of gaussian smoothing to the RT and DT dimensions.
        factors (list of Factor objects): List of Factor objects resulting from the factorize method.

    """

    def __init__(self, source_file, tensor_idx, timepoint_idx, name,
                 total_mass_window, n_concatenated, charge_states, integrated_mz_limits, bins_per_isotope_peak, normalization_factor, **kwargs):
        """Initializes an instance of the DataTensor class from 

        Args:
            source_file (str): Path of DataTensor's parent resources/tensors/.cpickle.zlib file.
            tensor_idx (int): Deprecated - Index of this tensor in a concatenated tensor, now always 0.
            timepoint_idx (int): Index of tensor's source timepoint in config["timepoints"].
            name (str): Name of rt-group DataTensor is a member of.
            total_mass_window (int): Magnitude of DataTensor's m/Z dimension in bins.
            n_concatenated (int): Number of DataTensors combined to make the instance DataTensor. Deprecated value always 1. 
            charge_states (list with one int): Net positive charge on protein represented in DataTensor, format left over from concatenation.
            integrated_mz_limits (numpy 2D array): Values for low and high m/Z limits of integration around expected peak centers.
            bins_per_isotope_peak (int):  Number of integrated m/Z bins for a signal to be considered as an IsotopeCluster.
            normalization_factor (float): Divisor for the integrated m/Z intensity of any IsotopeCluster from a parent DataTensor instance. 

        Keyword Arguments: 
            rts (numpy array): Intensity summed over each rt bin, from source tensor file.
            dts (numpy array): Intensity summer over each dt bin, from source tensor file.
            seq_out (numpy array): Intensity summed over each m/Z bin in a flat array, from source tensor file.
            int_seq_out (numpy array): Intensities integrated over all dimensions in a flat array. (TODO: Review, not sure about this)
            concat_dt_idxs (list of ints): Deprecated - Indices marking the boundaries between different DataTensors' concatenated DT dimensions.
            concatenated_grid (numpy array): Deprecated - int_grid_out from multiple DataTensors concatenated along the dt axis.

        """
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

        # Normal case, single tensor.
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

        # Handle concatenated tensor case, check for required inputs TODO: Remove? Deprecated, but Gabe has said to keep this before, ask again.
        else:
            if not all("dts" in kwargs and "rts" in kwargs and
                       "lows" in kwargs and "highs" in kwargs and
                       "concatenated_grid" in kwargs and
                       "abs_mz_low" in kwargs and "concat_dt_idxs" in kwargs):

                print("Concatenated Tensor Missing Required Values")
                sys.exit()


    def sparse_to_full_tensor_reprofile(self, data, integrated_mz_limits, bins_per_isotope_peak = 7, ms_resolution=25000):
        """Takes a tuple of three vectors and performs the 3-vector outer-product and applies resolution scaling.

        Args:
            data (tuple of 3 list-likes): Tuple of 3 vectors whose outer product recreates the extracted signal.
            integrated_mz_limits (list of list-likes of floats): The upper and lower bounds on integration surrounding each
                m/Z-peak center, a list of list-likes of length 2.
            bins_per_isotope_peak (int): Minimum number of bins a signal must span to be considered as an IsotopeCluster.
            ms_resolution (int): Number of bins in the m/Z dimension of the reprofiled tensor.

        Returns:
            retention_labels (list of floats): List of the absolute LC-RT values associated with each RT bin.
            drift_labels (list of floats):  List of the absolute IMS-DT values associated with each DT bin.
            mz_bin_centers (list of floats):  List of the absolute m/Z values associated with each m/Z bin.
            bins_per_isotope_peak (int): Minimum number of bins a signal must span to be considered as an IsotopeCluster.
            tensor3_out (numpy array): Full and reprofiled 3D tensor of LC-IMS-MS data.

        """
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


    def gauss(self, grid, rt_sig=3, dt_sig=1):
        """Applies a gaussian filter to the first two dimensions of a 3D tensor.

        Args:
            grid (numpy array): 3D tensor of LC-IMS-MS data, self.full_grid_out.
            rt_sig (int or float): Gaussian sigma for smoothing function on LC-RT dimension, default 3.
            dt_sig (int or float): Gaussian sigma for smoothing function on IMS-DT dimension, defualt 1.

        Returns:
            gauss_grid (type): Input tensor after applying gaussian smoothing.

        """

        gauss_grid = np.zeros(np.shape(grid))
        for i in range(np.shape(grid)[2]):
            gauss_grid[:, :, i] = gaussian_filter(grid[:, :, i],
                                                  (rt_sig, dt_sig))
        return gauss_grid


    # TODO: This isn't great style, make this take the tensor as input and return the factors.
    def factorize(self, n_factors=4, new_mz_len=None, gauss_params=None): 
        """Performs the non-negative PARAFAC on tensor, implemented by nn-fac python library, saves to self.Factors.

        Args:
            n_factors (int): The number of factors used to decompose the input tensor.
            new_mz_len (int): Number of bins desired in output tensor m/Z dimension, performs interpolation.
            gauss_params (tuple of 2 ints): Two values indicating the width of smoothing in LC-RT and IMS-DT dimensions respectively.

        Returns:
            None
        """
        # Test factorization starting at n_factors = 15 and counting down, keep factorization that has no factors with correlation greater than 0.2 in any dimension.

        def corr_check(factors):
            """Checks the maximum correlation between Factors in each dimension, used to determine if n_factors should be reduced.

            Args:
                factors (list of Factor objects): Factors resulting from PARAFAC being checked for inter-correlation.

            Returns:
                maximum_correlation (float): Maximum correlation between any two factors in any dimension.

            """
            # Checks scipy non_negatve_parafac output factors for inter-factor (off-diagonal) correlations > cutoff, returns True if all values are < cutoff

            a = np.minimum(
                np.minimum(np.corrcoef(factors[0].T),
                           np.corrcoef(factors[1].T)),
                np.corrcoef(factors[2].T),
            )

            return np.max(a[np.where(~np.eye(a.shape[0], dtype=bool))])


        def pmem(id_str):
            """Prints memory in use by process with a passed debug label.

            Args:
                id_str (str): String to prepend to memory output, for identifying origin of pmem call.

            Returns:
                None

            """
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
    """Computes the cumulative distribution function for a gaussian with given center, bounds, and width.

    Args:
        low_bound (float): Lower or left bound on gaussian, unitless. 
        upper_bound (float):  Upper or right bound on gaussian, unitless.
        center (float): Center of gaussian, unitless. 
        width (float): Width of gaussian, unitless.

    Returns:
        auc (float): Area Under the Curve, computed cumulative distribution of specified gaussian.

    """
    lb_cdf = norm.cdf(low_bound, loc=center, scale=width)
    ub_cdf = norm.cdf(upper_bound, loc=center, scale=width)
    auc = ub_cdf - lb_cdf
    return auc


def estimate_gauss_param(y_data, x_data):
    """Estimates the parameters of a gaussian fit to a set of datapoints with x,y coordinates in x_data and y_data.

    Args:
        y_data (iterable of floats): Y values of data to fit.
        x_data (iterable of floats): X values of data to fit.

    Returns:
        init_guess (list of floats): baseline offset, amplitude, center, width.

    """
    ymax = np.max(y_data)
    maxindex = np.nonzero(y_data == ymax)[0]
    peakmax_x = x_data[maxindex][0]
    norm_arr = y_data/max(y_data)
    bins_for_width = norm_arr[norm_arr > 0.70]
    width_bin = len(bins_for_width)
    init_guess = [0, ymax, peakmax_x, width_bin]
    # bounds = ([0, 0, 0, 0], [np.inf, np.inf, len(x_data)-1, len(x_data)-1])
    return init_guess


def gauss_func(x, y0, A, xc, w):
    """Model Gaussian function to pass to scipy.optimize.curve_fit.

    Args:
        x (list of floats): X dimension values.
        y0 (float): Offset of y-values from 0.
        A (float): Amplitude of Gaussian.
        xc (float): Center of Gaussian in x dimension.
        w (float): Width or sigma of Gaussian. 

    Returns:
        y (list of floats): Y-values of Gaussian function evaluated over x.

    """
    rxc = ((x - xc) ** 2) / (2 * (w ** 2))
    y = y0 + A * (np.exp(-rxc))
    return y


def adjrsquared(r2, param, num):
    """Calculates the adjusted R^2 for a parametric fit to data.

    Args:
         r2 (float): R^2 or 'coefficient of determination' of a linear regression between fitted and observed values.
         param (int): Number of parameters in fitting model.
         num (int): Number of datapoints in sample.

    Returns:
        y (float): Adjusted R^2 <= R^2, increases when additional parameters improve fit more than could be expected by chance.

    """
    y = 1 - (((1 - r2) * (num - 1)) / (num - param - 1))
    return y


def fit_gaussian(x_data, y_data, data_label="dt"):
    """Performs fitting of a gaussian function to a provided data sample and computes linear regression on residuals to measure quality of fit.

    Args:
        x_data (list of floats): X dimension values for sample data.
        y_data (list of floats): Y dimension values for sample data.
        data_label (str): Label indicating origin of data in multidimensional separation (rt, dt, m/z). 

    Returns:
        gauss_fit_dict (dict): Contains the following key-value pairs describing the Guassian and linear regression parameters.
            "gauss_fit_success" (bool): Boolean indicating the success (True) or failure (False) of the fitting operation.
            "y_baseline" (float): Fitted parameter for the Gaussian function's offset from y=0. 
            "y_amp" (float): Fitted parameter for the amplitude of the Gaussian function.
            "xc" (float): Fitted parameter for the center of the Gaussian in the x dimension.
            "width" (float): Fitted parameter for the x dimensional width of the Gaussian function.
            "y_fit" (list of floats): Y values of fitted Gaussian function evaluated over x. 
            "fit_rmse" (float): Root-mean-square error, the standard deviation of the residuals between the fit and sample.
            "fit_lingress_slope" (float): The slope of the linear regression line over the residuals between fit and sample.
            "fit_lingress_intercept" (float): The intercept point of the line fit to the residuals.
            "fit_lingress_pvalue" (float): The p-value for a hypothesis test whose null hypothesis is that the above slope is zero
            "fit_lingress_stderr" (float): Standard error of the estimated slope under the assumption of residual normality.
            "fit_lingress_r2" (float): R^2 or 'coeffiecient of determination' of linear regression over residuals.
            "fit_lingress_adj_r2" (float): Adjusted R^2, always <= R^2, decreases with extraneous parameters.
            "auc" (float): Area under the curve, cumulative distribution function of fitted gaussian evaluated over the length of x_data.

    """
    init_guess = estimate_gauss_param(y_data, x_data)
    gauss_fit_dict = dict()
    gauss_fit_dict['data_label'] = data_label
    gauss_fit_dict['gauss_fit_success'] = False
    gauss_fit_dict['xc'] = center_of_mass(y_data)[0]
    gauss_fit_dict['auc'] = 1.0
    gauss_fit_dict['fit_rmse'] = 100.0
    gauss_fit_dict['fit_linregress_r2'] = 0.0
    gauss_fit_dict['fit_lingress_adj_r2'] = 0.0

    try:
        popt, pcov = curve_fit(gauss_func, x_data, y_data, p0=init_guess, maxfev=100000)
        if popt[2] < 0:
            return gauss_fit_dict
        if popt[3] < 0:
            return gauss_fit_dict
        else:
            y_fit = gauss_func(x_data, *popt)
            fit_rmse = mean_squared_error(y_data/max(y_data), y_fit/max(y_fit), squared=False)
            slope, intercept, rvalue, pvalue, stderr = linregress(y_data, y_fit)
            adj_r2 = adjrsquared(r2=rvalue**2, param=4, num=len(y_data))
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
            gauss_fit_dict['auc'] = cal_area_under_curve_from_normal_distribution(low_bound=x_data[0],
                                                                                  upper_bound=x_data[-1],
                                                                                  center=popt[2],
                                                                                  width=popt[3])
            return gauss_fit_dict
    except:
        return gauss_fit_dict


def model_data_with_gauss(x_data, gauss_params):
    # TODO: Add docstring.

    data_length = len(x_data)
    bin_value = (x_data[-1] - x_data[0])/data_length
    center = gauss_params[2]
    data_length_half = int(data_length/2)
    low_val = center - (data_length_half * bin_value)
    new_x_data = []
    for num in range(data_length):
        val = low_val + (num * bin_value)
        new_x_data.append(val)
    new_x_data = np.array(new_x_data)
    new_gauss_data = gauss_func(new_x_data, *gauss_params)
    return new_gauss_data


class Factor:
    """A container for a factor output from DataTensor.factorize(), may contain the isolated signal from a single charged species.

    If the class has public attributes, they may be documented here
    in an ``Attributes`` section and follow the same formatting as a
    function's ``Args`` section. Alternatively, attributes may be documented
    inline with the attribute's declaration (see __init__ method below).

    Attributes:
        source_file (str): Path of DataTensor's parent resources/tensors/.cpickle.zlib file.
        tensor_idx (int): Index of Factor's parent DataTensor in a concatenated DataTensor. Deprecated, this value always 0.
        timepoint_idx (int): Index of the Factor's HDX timepoint in config["timepoints"].
        name (str): Name of Factor's rt-group.
        charge_states (list with one int): Net positive charge on protein represented in parent DataTensor.
        integrated_mz_limits (numpy 2D array): Values for low and high m/Z limits of integration around expected peak centers.
        rts (numpy array): Intensity of Factor summed over each rt bin.
        dts (numpy array): Intensity of Factor summed over each dt bin.
        mz_data (numpy array): Intensity of Factor summed over each m/Z bin in a flat array.
        retention_labels (list of floats): Mapping of RT bins to corresponding absolute retention time in minutes.
        drift_labels (list of floats): Mapping of DT bins to corresponding absolute drift time in miliseconds.
        mz_labels (list of floats): Mapping of m/Z bins to corresponding m/Z.
        factor_idx (int): Index of instance Factor in DataTensor.factors list.
        bins_per_isotope_peak (int): Number of integrated m/Z bins for a signal to be considered as an IsotopeCluster.
        n_concatenated (int): Number of DataTensors combined to make the instance DataTensor. Deprecated, value always 1.
        concat_dt_idxs (list of ints): Deprecated - Indices marking the boundaries between different DataTensors' concatenated DT dimensions.
        normalization_factor (float): Divisor for the integrated m/Z intensity of any IsotopeCluster from a parent DataTensor instance.
        integrated_mz_data (numpy array): Coarsened version of mz_data. Number of bins to sum per index is determined by bins_per_isotope_peak.
        max_rtdt (float): Product of maximal values from rts and dts.
        outer_rtdt (float): Sum of outer product of rts and dts, multiplied by the sum of mz_data to find the magnitude of the Factor.
        integrated_mz_baseline (numpy array): Baseline signal to subtract from mz_data. Deprecated.
        baseline_subtracted_integrated_mz (numpy array): Deprecated, copy of integrated_mz_data.
        rt_gauss_fit_success (bool): Boolean indicating success of Gaussian fit operation on rts. True => success, False => failure.
        rt_auc (float): Cumulative distribution function of Gaussian fit to rts evaluated between estimated bounds.
        rt_com (float): Computed center-of-mass of the Gaussian fit to rts. 
        rt_gaussian_rmse (float): Root-mean-square error, the standard deviation of the residuals between fitted values and rts.
        rt_gauss_fit_r2 (float): R^2 or 'coeffiecient of determination' of linear regression over residuals between fitted values and rts.
        dt_gauss_fit_success (bool): Boolean indicating success of Gaussian fit operation on dts. True => success, False => failure.
        dt_auc (float): Cumulative distribution function of Gaussian fit to dts evaluated between estimated bounds.
        dt_com (float): Computed center-of-mass of the Gaussian fit to dts. 
        dt_gaussian_rmse (float): Root-mean-square error, the standard deviation of the residuals between fitted values and dts.
        dt_gauss_fit_r2 (float): R^2 or 'coeffiecient of determination' of linear regression over residuals between fitted values and dts.
        isotope_clusters (list of IsotopeCluster objects): Contains IsotopeCluster objects made from candidate signals in integrated_mz_data

    """
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
        """Creates an instance of the Factor class from one factor of a PARAFAC run.

        Args:
            source_file (str): Path of DataTensor's parent resources/tensors/.cpickle.zlib file.
            tensor_idx (int): Index of Factor's parent DataTensor in a concatenated DataTensor. Deprecated, this value always 0.
            timepoint_idx (int): Index of the Factor's HDX timepoint in config["timepoints"].
            name (str): Name of Factor's rt-group.
            charge_states (list with one int): Net positive charge on protein represented in parent DataTensor.
            integrated_mz_limits (numpy 2D array): Values for low and high m/Z limits of integration around expected peak centers.
            rts (numpy array): Intensity of Factor summed over each rt bin.
            dts (numpy array): Intensity of Factor summed over each dt bin.
            mz_data (numpy array): Intensity of Factor summed over each m/Z bin in a flat array.
            retention_labels (list of floats): Mapping of RT bins to corresponding absolute retention time in minutes.
            drift_labels (list of floats): Mapping of DT bins to corresponding absolute drift time in miliseconds.
            mz_labels (list of floats): Mapping of m/Z bins to corresponding m/Z.
            factor_idx (int): Index of instance Factor in DataTensor.factors list.
            bins_per_isotope_peak (int): Number of integrated m/Z bins for a signal to be considered as an IsotopeCluster.
            n_concatenated (int): Number of DataTensors combined to make the instance DataTensor. Deprecated, value always 1.
            concat_dt_idxs (list of ints): Deprecated - Indices marking the boundaries between different DataTensors' concatenated DT dimensions.
            normalization_factor (float): Divisor for the integrated m/Z intensity of any IsotopeCluster from a parent DataTensor instance.

        """
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

        # Integrate within expected peak bounds.
        self.integrated_mz_data = np.sum(np.reshape(mz_data, (-1, self.bins_per_isotope_peak)), axis=1)

        # This can be a shared function
        #self.integrated_mz_baseline = peakutils.baseline(
        #    np.asarray(self.integrated_mz_data),
        #    6)  # 6 degree curve seems to work well
        #       
        #self.baseline_subtracted_integrated_mz = (self.integrated_mz_data -
        #                                          self.integrated_mz_baseline)
 
        # fit factor rts and dts to gaussian
        rt_gauss_fit = fit_gaussian(np.arange(len(self.rts)), self.rts, data_label='rt')
        dt_gauss_fit = fit_gaussian(np.arange(len(self.dts)), self.dts, data_label='dt')

        self.rt_gauss_fit_success = rt_gauss_fit['gauss_fit_success']
        self.rt_auc = rt_gauss_fit['auc']
        self.rt_com = rt_gauss_fit['xc']
        self.rt_gaussian_rmse = rt_gauss_fit['fit_rmse']
        self.rt_gauss_fit_r2 = rt_gauss_fit['fit_linregress_r2']

        self.dt_gauss_fit_success = dt_gauss_fit['gauss_fit_success']
        self.dt_auc = dt_gauss_fit['auc']
        self.dt_com = dt_gauss_fit['xc']
        self.dt_gaussian_rmse = dt_gauss_fit['fit_rmse']
        self.dt_gauss_fit_r2 = dt_gauss_fit['fit_linregress_r2']

        # calculate max rtdt and outer rtdt based on gauss fits
        if rt_gauss_fit['gauss_fit_success']:
            gauss_params = [rt_gauss_fit['y_baseline'], rt_gauss_fit['y_amp'], rt_gauss_fit['xc'], rt_gauss_fit['width']]
            rt_fac = model_data_with_gauss(np.arange(len(self.rts)), gauss_params)
        else:
            rt_fac = self.rts

        if dt_gauss_fit['gauss_fit_success']:
            gauss_params = [dt_gauss_fit['y_baseline'], dt_gauss_fit['y_amp'], dt_gauss_fit['xc'], dt_gauss_fit['width']]
            dt_fac = model_data_with_gauss(np.arange(len(self.dts)), gauss_params)
        else:
            dt_fac = self.rts

        # self.max_rtdt = max(self.rts) * max(self.dts)
        # self.outer_rtdt = sum(sum(np.outer(self.rts, self.dts)))

        self.max_rtdt_old = max(self.rts) * max(self.dts)
        self.outer_rtdt_old = sum(sum(np.outer(self.rts, self.dts)))

        self.max_rtdt = max(rt_fac) * max(dt_fac)
        self.outer_rtdt = sum(sum(np.outer(rt_fac, dt_fac)))

        # assign mean rt and dt values
        self.rt_mean = np.mean(np.arange(len(self.rts)))
        self.dt_mean = np.mean(np.arange(len(self.dts)))

        ## old protocol.
        # Writes to self.isotope_clusters
        # self.find_isotope_clusters()  # heuristic height value, should be high-level param TODO - Will require passage through DataTensor class

        ## now generates self.isotope_clusters upon calling the function self.find_isotope_clusters

    def rel_height_peak_bounds(self, centers, norm_integrated_mz, baseline_threshold=0.15, rel_ht_threshold=0.2):
        """Determines upper and lower bounds of IsotopeClusters within a Factor.

        Args:
            centers (type): List of peak positions in integrated m/Z from scipy.signal.find_peaks().
            norm_integrated_mz (numpy array): Integrated m/Z of Factor normalized for signal centeredness and timepoint variations.
            baseline_threshold (float): Minimum ratio of peak height to maximum height for a peak to be considered as an IC. 
            rel_ht_threshold (float): Ratio of point height to peak height that a point must exceed to be part of an IC.

        Returns:
            out (list of tuple pairs of ints): List of tuples with low and high bounding indices for ICs in the centers list.
        """
        out = []
        for center in centers:
            #REVIEW: Is this what this is supposed to be doing? Should be > max(norm_integrated_mz) * baseline_threshold? 
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

    def find_isotope_clusters(self, prominence=0.15, width_val=3, rel_height_filter=True, baseline_threshold=0.15, rel_height_threshold=0.10):
        """Identifies portions of the integrated mz dimension that look 'isotope-cluster-like', saves in isotope_clusters.

        Args:
            prominence (float): Ratio of array's maximum intesity that a peak must surpass to be considered.
            width_val (int): Minimum width to consider a peak to be an isotope cluster.
            rel_height_filter (bool): Switch to apply relative height filtering, True applies filter. 
            baseline_threshold (float): Minimum height of peak to consider it to be an isotope cluster.
            rel_height_threshold (float): Proportion determining the minimum intensity for bins near a peak to be part of an isotope cluster.

        Returns:
            None

        """
        self.isotope_clusters = []
        norm_integrated_mz = self.integrated_mz_data/max(self.integrated_mz_data)
        peaks, feature_dict = find_peaks(norm_integrated_mz,
                                         prominence=prominence,
                                         width=width_val)
        if len(peaks) == 0:
            ic_idxs = [(0, len(self.integrated_mz_data)-1)]
            int_mz_width = [2]
            #return
        else:
            int_mz_width = [
                feature_dict["widths"][i]
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

        cluster_idx = 0
        for integrated_indices, integrated_mz_width in zip(ic_idxs, int_mz_width):
            if integrated_indices != None:
                newIC = IsotopeCluster(
                    integrated_mz_peak_width=integrated_mz_width,
                    charge_states=self.charge_states,
                    factor_mz_data=copy.deepcopy(self.mz_data),
                    name=self.name,
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
                    rt_mean=self.rt_mean,
                    dt_mean=self.dt_mean,
                    rt_gauss_fit_success=self.rt_gauss_fit_success,
                    dt_gauss_fit_success=self.dt_gauss_fit_success,
                    rt_gaussian_rmse=self.rt_gaussian_rmse,
                    dt_gaussian_rmse=self.dt_gaussian_rmse,
                    rt_com=self.rt_com,
                    dt_coms=self.dt_com,
                    rt_auc=self.rt_auc,
                    dt_auc=self.dt_auc,
                    retention_labels=self.retention_labels,
                    drift_labels=self.drift_labels,
                    mz_labels=self.mz_labels,
                    bins_per_isotope_peak=self.bins_per_isotope_peak,
                    max_rtdt=self.max_rtdt,
                    outer_rtdt=self.outer_rtdt,
                    max_rtdt_old=self.max_rtdt_old,
                    outer_rtdt_old=self.outer_rtdt_old,
                    n_concatenated=self.n_concatenated,
                    concat_dt_idxs=self.concat_dt_idxs,
                    normalization_factor=self.normalization_factor
                )
                if (newIC.baseline_peak_error / newIC.baseline_auc <
                        0.2):  # TODO: HARDCODE
                    self.isotope_clusters.append(newIC)
                    cluster_idx += 1
        return


class IsotopeCluster:
    """Contains a portion of Factor.integrated_mz_data identified to have isotope-cluster-like characteristics, stores data of parent factor.

    Attributes:
        integrated_mz_peak_width (int): Number of isotope peaks estimated to be included in the isotope cluster. 
        charge_states (list with one int): Net positive charge on protein represented in parent DataTensor.
        factor_mz_data (numpy array): The full mz_data of the parent Factor.
        name (str): Name of IC's rt-group.
        source_file (str): Path of DataTensor's parent resources/tensors/.cpickle.zlib file.
        tensor_idx (int): Index of the IC's parent DataTensor in a concatenated DataTensor. Deprecated, this value always 0.
        timepoint_idx (int): Index of the IC's HDX timepoint in config["timepoints"].
        n_factors (int): Number of factors used in decomposition of IC's parent DataTensor.
        factor_idx (int): Index of IC's parent Factor in its parent DataTensor.factors.
        cluster_idx (int): Index of IC in parent Factor's Factor.isotope_clusters.
        low_idx (int): Lower bound index of IC in integrated m/Z dimension.
        high_idx (int): Upper bound index of IC in integrated m/Z dimension.
        rts (numpy array): Intensity of IC's parent Factor summed over each rt bin.
        dts (numpy array): Intensity of IC's parent Factor summed over each dt bin.
        rt_mean (float):
        dt_mean (float):
        rt_gauss_fit_success (bool): Boolean indicating success of Gaussian fit operation on rts. True => success, False => failure.
        rt_auc (float): Cumulative distribution function of Gaussian fit to rts evaluated between estimated bounds.
        rt_com (float): Computed center-of-mass of the Gaussian fit to rts. 
        rt_gaussian_rmse (float): Root-mean-square error, the standard deviation of the residuals between fitted values and rts.
        rt_gauss_fit_r2 (float): R^2 or 'coeffiecient of determination' of linear regression over residuals between fitted values and rts.
        dt_gauss_fit_success (bool): Boolean indicating success of Gaussian fit operation on dts. True => success, False => failure.
        dt_auc (float): Cumulative distribution function of Gaussian fit to dts evaluated between estimated bounds.
        dt_com (float): Computed center-of-mass of the Gaussian fit to dts. 
        dt_gaussian_rmse (float): Root-mean-square error, the standard deviation of the residuals between fitted values and dts.
        dt_gauss_fit_r2 (float): R^2 or 'coeffiecient of determination' of linear regression over residuals between fitted values and dts.
        retention_labels (list of floats): Mapping of RT bins to corresponding absolute retention time in minutes.
        drift_labels (list of floats): Mapping of DT bins to corresponding absolute drift time in miliseconds.
        mz_labels (list of floats): Mapping of m/Z bins to corresponding m/Z.
        bins_per_isotope_peak (int): Number of integrated m/Z bins for a signal to be considered as an IsotopeCluster.
        max_rtdt (float):
        outer_rt_dt (float):
        max_rtdt_old(float):
        outer_rtdt_old (float):
        n_concatenated (int): Number of DataTensors combined to make the instance DataTensor. Deprecated, value always 1.
        concat_dt_idxs (list of ints): Deprecated - Indices marking the boundaries between different DataTensors' concatenated DT dimensions.
        normalization_factor (float): Divisor for the integrated m/Z intensity of any IsotopeCluster from a parent DataTensor instance.
        cluster_mz_data (numpy array): 
        auc (float): 
        baseline (float): Deprecated, always 0
        baseline_subtracted_mz (numpy array): Deprecated, always a copy of self.cluster_mz_data.
        baseline_auc (float): Deprecated, always a copy of auc.
        log_baseline_auc (float): Base 10 log of baseline_auc.
        baseline_max_peak_height (float): Deprecated, max value of baseline_subtracted_mz

        # TODO: Document this line.
        peak_error (float): Average of the absolute distances of peak centers from expected centers of integration bounds.
        baseline_peak_error (float): Deprecated, copy of peak_error.

    """
    def __init__(
        self,
        integrated_mz_peak_width,
        charge_states,
        factor_mz_data,
        name,
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
        rt_mean,
        dt_mean,
        rt_gauss_fit_success,
        dt_gauss_fit_success,
        rt_gaussian_rmse,
        dt_gaussian_rmse,
        rt_com,
        dt_coms,
        rt_auc,
        dt_auc,
        retention_labels,
        drift_labels,
        mz_labels,
        bins_per_isotope_peak,
        max_rtdt,
        outer_rtdt,
        max_rtdt_old,
        outer_rtdt_old,
        n_concatenated,
        concat_dt_idxs,
        normalization_factor
    ):
        """Creates an instance of the IsotopeCluster class from a portion of Factor.mz_data.

        Args:
            integrated_mz_peak_width (int): Number of isotope peaks estimated to be included in the isotope cluster. 
            charge_states (list with one int): Net positive charge on protein represented in parent DataTensor.
            factor_mz_data (numpy array): The full mz_data of the parent Factor.
            name (str): Name of IC's rt-group.
            source_file (str): Path of DataTensor's parent resources/tensors/.cpickle.zlib file.
            tensor_idx (int): Index of the IC's parent DataTensor in a concatenated DataTensor. Deprecated, this value always 0.
            timepoint_idx (int): Index of the IC's HDX timepoint in config["timepoints"].
            n_factors (int): Number of factors used in decomposition of IC's parent DataTensor.
            factor_idx (int): Index of IC's parent Factor in its parent DataTensor.factors.
            cluster_idx (int): Index of IC in parent Factor's Factor.isotope_clusters.
            low_idx (int): Lower bound index of IC in integrated m/Z dimension.
            high_idx (int): Upper bound index of IC in integrated m/Z dimension.
            rts (numpy array): Intensity of IC's parent Factor summed over each rt bin.
            dts (numpy array): Intensity of IC's parent Factor summed over each dt bin.
            rt_mean (float):
            dt_mean (float):
            rt_gauss_fit_success (bool): Boolean indicating success of Gaussian fit operation on rts. True => success, False => failure.
            rt_auc (float): Cumulative distribution function of Gaussian fit to rts evaluated between estimated bounds.
            rt_com (float): Computed center-of-mass of the Gaussian fit to rts. 
            rt_gaussian_rmse (float): Root-mean-square error, the standard deviation of the residuals between fitted values and rts.
            rt_gauss_fit_r2 (float): R^2 or 'coeffiecient of determination' of linear regression over residuals between fitted values and rts.
            dt_gauss_fit_success (bool): Boolean indicating success of Gaussian fit operation on dts. True => success, False => failure.
            dt_auc (float): Cumulative distribution function of Gaussian fit to dts evaluated between estimated bounds.
            dt_com (float): Computed center-of-mass of the Gaussian fit to dts. 
            dt_gaussian_rmse (float): Root-mean-square error, the standard deviation of the residuals between fitted values and dts.
            dt_gauss_fit_r2 (float): R^2 or 'coeffiecient of determination' of linear regression over residuals between fitted values and dts.
            retention_labels (list of floats): Mapping of RT bins to corresponding absolute retention time in minutes.
            drift_labels (list of floats): Mapping of DT bins to corresponding absolute drift time in miliseconds.
            mz_labels (list of floats): Mapping of m/Z bins to corresponding m/Z.
            bins_per_isotope_peak (int): Number of integrated m/Z bins for a signal to be considered as an IsotopeCluster.
            max_rtdt (float):
            outer_rt_dt (float):
            max_rtdt_old(float):
            outer_rtdt_old (float):
            n_concatenated (int): Number of DataTensors combined to make the instance DataTensor. Deprecated, value always 1.
            concat_dt_idxs (list of ints): Deprecated - Indices marking the boundaries between different DataTensors' concatenated DT dimensions.
            normalization_factor (float): Divisor for the integrated m/Z intensity of any IsotopeCluster from a parent DataTensor instance.
        """
        self.integrated_mz_peak_width = integrated_mz_peak_width
        self.charge_states = charge_states
        self.factor_mz_data = factor_mz_data
        self.name = name
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
        self.rt_mean = rt_mean
        self.dt_mean = dt_mean
        self.rt_gauss_fit_success = rt_gauss_fit_success
        self.dt_gauss_fit_success = dt_gauss_fit_success
        self.rt_gaussian_rmse = rt_gaussian_rmse
        self.dt_gaussian_rmse = dt_gaussian_rmse
        self.rt_com = rt_com
        self.dt_coms = dt_coms
        self.rt_auc = rt_auc
        self.dt_auc = dt_auc
        self.retention_labels = retention_labels
        self.drift_labels = drift_labels
        self.mz_labels = mz_labels
        self.bins_per_isotope_peak = bins_per_isotope_peak
        self.max_rtdt_old = max_rtdt_old
        self.outer_rtdt_old = outer_rtdt_old
        self.max_rtdt = max_rtdt
        self.outer_rtdt = outer_rtdt
        self.n_concatenated = n_concatenated
        self.concat_dt_idxs = concat_dt_idxs
        self.normalization_factor = normalization_factor

        # Prune factor_mz to get window around cluster that is consistent between charge-states.
        self.cluster_mz_data = copy.deepcopy(self.factor_mz_data)
        self.cluster_mz_data[0:self.low_idx] = 0
        self.cluster_mz_data[self.high_idx:] = 0

        # Integrate area of IC and normalize according the TIC counts.
        self.auc = sum(self.cluster_mz_data) * self.outer_rtdt / self.normalization_factor

        # Reshape cluster m/Z data by expected number of bins per isotope.
        isotope_peak_array = np.reshape(self.cluster_mz_data, (-1, self.bins_per_isotope_peak))
        self.baseline = 0
        self.baseline_subtracted_mz = self.cluster_mz_data 
        self.baseline_auc = self.auc
        self.log_baseline_auc = np.log(self.baseline_auc)
        self.baseline_max_peak_height = max(self.baseline_subtracted_mz)
        self.baseline_integrated_mz = np.sum(isotope_peak_array, axis=1)

        # Takes the average of the absolute distances of peak centers from expected centers of integration bounds.
        self.peak_error = np.average(np.abs(np.argmax(isotope_peak_array,axis=1) - ((self.bins_per_isotope_peak - 1)/2)) 
                            / ((self.bins_per_isotope_peak - 1)/2), weights=self.baseline_integrated_mz)
        self.baseline_peak_error = self.peak_error

        # Cache int_mz and rt scoring values.
        # Compute baseline_integrated_mz_com, baseline_integrated_mz_std, baseline_integrated_mz_FWHM, and baseline_integrated_mz_rmse from gaussian fit.
        # Define functions for Gaussian fit.
        def gaussian_function(x, H, A, x0, sigma):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


        def gauss_fit(x, y):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            mean = sum(x * y) / sum(y)
            sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
            nonzeros = [index for index, value in enumerate(list(y)) if value!=0]
            popt, pcov = curve_fit(gaussian_function, x, y, p0=[0, max(y), mean, sigma], bounds=([0, 0, nonzeros[0], 0], [np.inf, np.inf, nonzeros[-1], np.inf]))
            return popt


        def params_from_gaussian_fit(self):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            try:
                x_data = [i for i in range(len(self.baseline_integrated_mz))]
                y_data = self.baseline_integrated_mz
                H, A, x0, sigma = gauss_fit(x_data, y_data)
                y_gaussian_fit = gaussian_function(x_data, *gauss_fit(x_data, y_data))
                rmse = mean_squared_error(y_data/max(y_data), y_gaussian_fit/max(y_gaussian_fit), squared=False)
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

        # rt_com is pre calculated during factor data class
        # self.rt_com = center_of_mass(self.rts)[0]

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

            ## dt coms is pre calculated during factor data class
            # self.dt_coms = [center_of_mass(dt)[0] for dt in single_dts]

            self.dt_norms = [dt / np.linalg.norm(dt) for dt in single_dts]
        else:
            ## dt coms is pre calculated during factor data class
            # self.dt_coms = [center_of_mass(self.dts)[0]]
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
