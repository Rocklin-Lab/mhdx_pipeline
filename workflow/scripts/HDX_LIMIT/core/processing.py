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
import sys
import psutil
import copy
import math
import molmass
import numpy as np
import pandas as pd
from scipy.stats import gmean, norm, linregress
from HDX_LIMIT import io, datatypes
from numpy import linspace, cumsum, searchsorted
from HDX_LIMIT.plot_factor_data import plot_factor_data_from_data_dict, plot_factor_data_from_data_tensor

from bokeh.plotting import figure
from bokeh.palettes import Spectral6
from bokeh.transform import linear_cmap
from bokeh.layouts import gridplot, column
from bokeh.models import HoverTool, ColorBar, Text, Div, Whisker
from bokeh.models.glyphs import MultiLine, Line
from bokeh.io import save, output_file
from bokeh.models.callbacks import CustomJS
from bokeh.models.sources import ColumnDataSource, CDSView
from bokeh.models.filters import Filter, GroupFilter, IndexFilter

import scipy as sp
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error


def filter_factors_on_rt_dt_gauss_fit(factor_list, rt_r2_cutoff=0.90, dt_r2_cutoff=0.90):
    """Filter Factors based on quality of RT and DT gaussian fit.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

    """

    """
    New factor list is created if the Factor 
    DT and RT gaussian fit r^2 values are high. If none of the Factors pass
    the filtering criteria, returns original Factor list.
    :param factor_list: gauss fitted factor list
    :param rt_r2_cutoff: rt gauss fit r2 cutoff
    :param dt_r2_cutoff: dt gauss fit r2 cutoff
    :return: filtered factor list
    """

    filtered_factors = []

    for factor in factor_list:
        if factor.rt_gauss_fit_r2 >= rt_r2_cutoff:
            if factor.dt_gauss_fit_r2 >= dt_r2_cutoff:
                filtered_factors.append(factor)

    new_factor_list = filtered_factors
    if len(filtered_factors) == 0:
        new_factor_list = factor_list

    return new_factor_list


def create_factor_data_object(data_tensor, gauss_params, timepoint_label=None):
    """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

    """

    """
    function to store factor data to factor data class
    :param data_tensor:
    :param gauss_params:
    :param timepoint_label:
    :return:
    """

    factor_data_dict = {
    "name": data_tensor.DataTensor.name,
    "charge_state": data_tensor.DataTensor.charge_states[0],
    "timepoint_index": data_tensor.DataTensor.timepoint_idx,
    "timepoint_label": timepoint_label,
    "retention_labels": data_tensor.DataTensor.retention_labels,
    "drift_labels": data_tensor.DataTensor.drift_labels,
    "mz_labels": data_tensor.DataTensor.mz_labels,
    "bins_per_isotope_peak": data_tensor.DataTensor.bins_per_isotope_peak,
    "tensor_3d_grid": data_tensor.DataTensor.full_grid_out,
    "gauss_params": gauss_params,
    "num_factors": len(data_tensor.DataTensor.factors),
    "factors": []
    }

    for num, factor in enumerate(data_tensor.DataTensor.factors):
        factor_dict = {
        "factor_num": factor.factor_idx,
        "factor_dt": factor.dts,
        "factor_rt": factor.rts,
        "factor_mz": factor.mz_data,
        "factor_integrated_mz": factor.integrated_mz_data
        }
        factor_data_dict['factors'].append(factor_dict)

    return factor_data_dict


def generate_tensor_factors(tensor_fpath, library_info_df, timepoint_index, gauss_params, mz_centers, normalization_factors,
                            n_factors=15,
                            factor_output_fpath=None,
                            factor_plot_output_path=None,
                            timepoint_label=None,
                            filter_factors=False,
                            factor_rt_r2_cutoff=0.90,
                            factor_dt_r2_cutoff=0.90):
    """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

    """

    """
    generate data tensor from a given tensor file path, library info file path, and timepoint index
    :param tensor_fpath: tensor file path
    :param library_info_fpath: library info file path
    :param timepoint_index: timepoint index int
    :param gauss_params: gaussian filter params in tuple (rt_sigma, dt_sigma)
    :param n_factors: maximum number of factors to be considered.
    :param timepoint_label: timepoint label (in sec, hr, min, etc)
    :return: data_tensor
    """

    #memory calculations
    process = psutil.Process(os.getpid())
    # memory before init
    print("Pre-Tensor-Initialization: " + str(process.memory_info().rss /
                                       (1024 * 1024 * 1024)))

    # data tensor initialization
    data_tensor = TensorGenerator(filename=tensor_fpath,
                                  library_info=library_info_df,
                                  timepoint_index=timepoint_index,
                                  mz_centers=mz_centers,
                                  normalization_factors=normalization_factors)

    print("Post-Tensor-Pre-Factor-Initialization: " + str(process.memory_info().rss /
                                        (1024 * 1024 * 1024)))

    print('Factorizing ... ')

    data_tensor.DataTensor.factorize(n_factors=n_factors,
                                     gauss_params=gauss_params)

    # profile memory after factorization
    print("Post-Factorization: " + str(process.memory_info().rss /
                                       (1024 * 1024 * 1024)))

    if filter_factors:
        filtered_factors = filter_factors_on_rt_dt_gauss_fit(factor_list=data_tensor.DataTensor.factors,
                                                             rt_r2_cutoff=factor_rt_r2_cutoff,
                                                             dt_r2_cutoff=factor_dt_r2_cutoff)
        data_tensor.DataTensor.factors = filtered_factors

    # save factor data object
    if factor_output_fpath != None:
        # create factor data dictionary
        factor_data_dictionary = create_factor_data_object(data_tensor=data_tensor,
                                                           gauss_params=gauss_params,
                                                           timepoint_label=timepoint_label)
        io.limit_write(factor_data_dictionary, factor_output_fpath)

    # plot_factor_data
    if factor_plot_output_path != None:
        plot_factor_data_from_data_tensor(data_tensor=data_tensor,
                                          output_path=factor_plot_output_path)

    return data_tensor


class TensorGenerator:
    """The summary line for a class docstring should fit on one line.

    If the class has public attributes, they may be documented here
    in an ``Attributes`` section and follow the same formatting as a
    function's ``Args`` section. Alternatively, attributes may be documented
    inline with the attribute's declaration (see __init__ method below).

    Properties created with the ``@property`` decorator should be documented
    in the property's getter method.

    Attributes:
        attr1 (str): Description of `attr1`.
        attr2 (:obj:`int`, optional): Description of `attr2`.

    """
    hd_mass_diff = 1.006277
    c13_mass_diff = 1.00335


    def __init__(self, filename, timepoint_index, library_info, mz_centers, normalization_factors, **kwargs):
        """Example of docstring on the __init__ method.

        The __init__ method may be documented in either the class level
        docstring, or as a docstring on the __init__ method itself.

        Either form is acceptable, but the two should not be mixed. Choose one
        convention to document the __init__ method and be consistent with it.

        Note:
            Do not include the `self` parameter in the ``Args`` section.

        Args:
            param1 (str): Description of `param1`.
            param2 (:obj:`int`, optional): Description of `param2`. Multiple
                lines are supported.
            param3 (:obj:`list` of :obj:`str`): Description of `param3`.

        """

        ###Set Instance Attributes###

        self.filename = filename
        self.timepoint_index = timepoint_index
        self.library_info = library_info
        self.mz_centers = mz_centers
        self.normalization_factors = normalization_factors
        my_mzml = ".".join("_".join(self.filename.split("_")[-5:]).split(".")[:2]) # Fix for updated rt-group directories.
        print(my_mzml)
        self.normalization_factor = normalization_factors.loc[normalization_factors["mzml"]==my_mzml]["normalization_factor"].values[0]

        if (
                kwargs is not None
        ):  
            for key in kwargs.keys():
                setattr(self, key, kwargs[key])

        if not hasattr(self, "low_mass_margin"):
            self.low_mass_margin = 10
        if not hasattr(self, "high_mass_margin"):
            self.high_mass_margin = 17
        if not hasattr(self, "ppm_radius"):
            self.ppm_radius = 30
        if not hasattr(self, "n_factors_low"):
            self.n_factors_low = 1
        if not hasattr(self, "n_factors_high"):
            self.n_factors_high = 3
        if not hasattr(self, "gauss_params"):
            self.gauss_params = (3, 1)
        if not hasattr(self, "bins_per_isotope_peak"):
            self.bins_per_isotope_peak = 7

        self.tensor = io.limit_read(self.filename)
        self.name = filename.split("/")[-2] # Expects format: path/to/{rt-group-name}/{rt-group-name}_{charge}_{file.mzML.gz}.cpickle.zlib.
        self.charge = int(filename.split("_")[-6]) # This is bad but works for the naming scheme.
        self.lib_idx = self.library_info.loc[(library_info["name"]==self.name) & (library_info["charge"]==self.charge)].index
        self.max_peak_center = len(self.library_info.loc[
            self.library_info["name"] == self.name]["sequence"].values[0])
        self.total_isotopes = self.max_peak_center + self.high_mass_margin
        self.total_mass_window = self.low_mass_margin + self.total_isotopes

        i = self.lib_idx
        self.mz_lows = self.library_info["obs_mz"].values[i] - (
            self.low_mass_margin / self.library_info["charge"].values[i])
        self.mz_highs = self.library_info["obs_mz"].values[i] + (
            self.total_isotopes / self.library_info["charge"].values[i])

        low_mz_limits = [center * ((1000000.0 - self.ppm_radius) / 1000000.0) for center in self.mz_centers]
        high_mz_limits = [center * ((1000000.0 + self.ppm_radius) / 1000000.0) for center in self.mz_centers]

        self.integrated_mz_limits = np.stack((low_mz_limits, high_mz_limits)).T

        # Instantitate DataTensor
        self.DataTensor = datatypes.DataTensor(
            source_file=self.filename,
            tensor_idx=self.lib_idx,
            timepoint_idx=self.timepoint_index,
            name=self.name,
            total_mass_window=self.total_mass_window,
            n_concatenated=1,
            charge_states=[self.library_info["charge"].values[self.lib_idx]],
            rts=self.tensor[0],
            dts=self.tensor[1],
            seq_out=self.tensor[2],
            int_seq_out=None,
            integrated_mz_limits=self.integrated_mz_limits,
            bins_per_isotope_peak=self.bins_per_isotope_peak,
            normalization_factor=self.normalization_factor
        )

        #self.DataTensor.lows = searchsorted(self.DataTensor.mz_labels,
        #                                    self.low_lims)
        #self.DataTensor.highs = searchsorted(self.DataTensor.mz_labels,
        #                                     self.high_lims)
        
        # Consider separating factorize from init
        # self.DataTensor.factorize(gauss_params=(3,1))


###
### Class - PathOptimizer:
### Optimizes HDX timeseries of Isotopic Clusters by generating a set of starting timeseries based on possible trajectory though the integrated mz dimension,
### Then iteratively uses a set of timeseries scoring functions to make the best single substitution until each starting series is optimized for score.
### Timeseries with best score at the end of all minimizations is selected as the winning path, which is output along with the alternatives for each timepoint.
###
class PathOptimizer:
    """The summary line for a class docstring should fit on one line.

    If the class has public attributes, they may be documented here
    in an ``Attributes`` section and follow the same formatting as a
    function's ``Args`` section. Alternatively, attributes may be documented
    inline with the attribute's declaration (see __init__ method below).

    Properties created with the ``@property`` decorator should be documented
    in the property's getter method.

    Attributes:
        attr1 (str): Description of `attr1`.
        attr2 (:obj:`int`, optional): Description of `attr2`.

    """

    """
    Generates sample 'paths' - trajectories through HDX timeseries - optimizes 'path' through hdx timeseries for all sample paths and selects an overall winning path.

    all_tp_clusters = <list> of <lists> of <TA.isotope_cluster>s for each HDX timepoint,

    """

    # TODO: add info_tuple-like struct, eventually change IC, PO, and bokeh related scoring systems to use dicts? Dicts would make changing column names simpler.

    def __init__(self,
                 name,
                 all_tp_clusters,
                 library_info,
                 timepoints,
                 n_undeut_runs,
                 prefilter=0,
                 old_data_dir=None,
                 **kwargs):
        """Example of docstring on the __init__ method.

        The __init__ method may be documented in either the class level
        docstring, or as a docstring on the __init__ method itself.

        Either form is acceptable, but the two should not be mixed. Choose one
        convention to document the __init__ method and be consistent with it.

        Note:
            Do not include the `self` parameter in the ``Args`` section.

        Args:
            param1 (str): Description of `param1`.
            param2 (:obj:`int`, optional): Description of `param2`. Multiple
                lines are supported.
            param3 (:obj:`list` of :obj:`str`): Description of `param3`.

        """

        # Set score weights
        self.int_mz_std_rmse_weight = 1
        self.baseline_peak_error_weight = 100
        self.delta_mz_rate_backward_weight = 1.65  # was 2
        self.delta_mz_rate_forward_weight = 1.65  # was 2
        self.int_mz_rot_fit_weight = 5
        self.dt_ground_fit_weight = 25
        self.rt_ground_fit_weight = 5
        self.rt_ground_rmse_weight = 10
        self.dt_ground_rmse_weight = 10
        self.auc_ground_rmse_weight = 20
        self.rmses_sum_weight = 1 
        self.maxint_sum_weight = 1
        self.int_mz_FWHM_rmse_weight = 1
        self.nearest_neighbor_penalty_weight = 1

        self.name = name
        self.all_tp_clusters = all_tp_clusters
        self.library_info = library_info
        self.prefilter = prefilter
        self.timepoints = timepoints
        self.n_undeut_runs = n_undeut_runs
        self.max_peak_center = len(
            self.library_info.loc[self.library_info["name"] ==
                                  self.name]["sequence"].values[0]
        )  # TODO: Fix. Bad solution, keep an eye on this, could break if only one charge in RT-group? Maybe always list, double check

        self.old_data_dir = old_data_dir
        self.old_files = None

        self.rt_com_cv = None
        self.dt_com_cv = None
        self.rt_error_rmse = None
        self.dt_error_rmse = None

        self.gather_old_data()
        self.select_undeuterated()
        self.precalculate_fit_to_ground()
        self.prefiltered_ics = self.weak_pareto_dom_filter()
        #self.prefiltered_ics = self.all_tp_clusters
        self.generate_sample_paths()


    def weak_pareto_dom_filter(self):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        self.all_tp_clusters[0] = [self.undeut_grounds[charge] for charge in self.undeut_grounds]
        out = []
        for tp in self.all_tp_clusters:
            tp_buffer = []
            for ic1 in tp:
                ic1_int_mz_com = np.round(ic1.baseline_integrated_mz_com)
                compare_flag = False
                for ic2 in tp:
                    if (
                        np.round(ic2.baseline_integrated_mz_com) == ic1_int_mz_com and
                        ic2.rt_ground_err**2 < ic1.rt_ground_err**2 and
                        ic2.dt_ground_err**2 < ic1.dt_ground_err**2 and
                        ic2.baseline_peak_error < ic1.baseline_peak_error and
                        ic2.rt_ground_fit > ic1.rt_ground_fit and 
                        ic2.dt_ground_fit > ic1.dt_ground_fit and 
                        ic2.auc_ground_err**2 < ic1.auc_ground_err**2
                       ):
                        compare_flag = True
                        break
                if not compare_flag:
                    tp_buffer.append(ic1)
            out.append(tp_buffer)
        return out


    def alt_weak_pareto_dom_filter(self):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        # Filters input of PO ICs to ensure no IC is worse in every score dimension than another IC (weak Pareto domination)

        out = []
        for tp in self.all_tp_clusters:

            tp_buffer = []
            center_dict = {}

            # make dict for tp int mz bins
            for i in range(len(tp[0].baseline_integrated_mz)):
                center_dict[i] = []

            # add ic to bin list closest to center
            for ic in tp:
                center_dict[min([len(tp[0].baseline_integrated_mz)-1,
                                   np.round(ic.baseline_integrated_mz_com)])].append(ic)

            # score all ics in each int_mz bin, keep only those that are not worse than another IC in all dimensions
            low_score_keys = ["rt_ground_err", "dt_ground_err", "peak_err"]
            high_score_keys = ["rt_ground_fit", "dt_ground_fit", "baseline_auc"]
            for i in range(len(tp[0].baseline_integrated_mz)):
                int_mz_buffer = []
                score_df = pd.DataFrame().from_dict({
                    "idx": [j for j in range(len(center_dict[i]))],
                    "rt_ground_err": [
                        abs(ic.rt_ground_err) for ic in center_dict[i]
                    ],
                    "dt_ground_err": [
                        abs(ic.dt_ground_err) for ic in center_dict[i]
                    ],
                    "peak_err": [
                        abs(ic.baseline_peak_error) for ic in center_dict[i]
                    ],
                    "rt_ground_fit": [
                        ic.rt_ground_fit for ic in center_dict[i]
                    ],
                    "dt_ground_fit": [
                        ic.dt_ground_fit for ic in center_dict[i]
                    ],
                    "baseline_auc": [ic.baseline_auc for ic in center_dict[i]],
                })

                if len(score_df) > 0:
                    for idx in score_df["idx"].values:
                        int_mz_dom_dict = dict.fromkeys(low_score_keys +
                                                        high_score_keys)

                        for key in low_score_keys:
                            dom_list = list(
                                score_df.sort_values(key)["idx"].values)
                            ic_pos = dom_list.index(idx)
                            int_mz_dom_dict[key] = set(dom_list[:ic_pos])

                        for key in high_score_keys:
                            dom_list = list(
                                score_df.sort_values(
                                    key, ascending=False)["idx"].values)
                            ic_pos = dom_list.index(idx)
                            int_mz_dom_dict[key] = set(dom_list[:ic_pos])

                        # Check if there is some IC with a better score in all dimensions by set intersection
                        if int_mz_dom_dict[low_score_keys[0]].intersection(
                                int_mz_dom_dict[low_score_keys[1]],
                                int_mz_dom_dict[low_score_keys[2]],
                                int_mz_dom_dict[high_score_keys[0]],
                                int_mz_dom_dict[high_score_keys[1]],
                        ):
                            # ic is weakly Pareto dominated, leave out of output
                            pass
                        else:
                            # ic is not weakly Pareto dominated, add to output
                            int_mz_buffer.append(center_dict[i][idx])

                for ic in int_mz_buffer:
                    tp_buffer.append(ic)

            out.append(tp_buffer)

        return out


    def gather_old_data(self):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        if self.old_data_dir is not None:
            self.old_files = sorted([
                fn for fn in glob.iglob(self.old_data_dir + "*.pickle")
                if "_".join(self.name.split("/")[-1].split("_")[:3]) in fn
            ])
            self.old_data = []
            for fn in self.old_files:
                ts = pickle.load(open(fn, "rb"))
                ts["charge"] = int(fn.split(".")[-3][-1])
                ts["delta_mz_rate"] = self.gabe_delta_mz_rate(ts["centroid"])
                ts["major_species_widths"] = [
                    len(np.nonzero(x)[0])
                    for x in ts["major_species_integrated_intensities"]
                ]
                self.old_data.append(ts)


    def select_undeuterated(self,
                            all_tp_clusters=None,
                            library_info=None,
                            name=None,
                            n_undeut_runs=None):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """

        """
        Selects undeuterated isotope cluster which best matches theoretically calculated isotope distribution for POI sequence, for each observed charge state of the POI
        all_tp_clusters = TensorGenerator attribute, e.g. T1 = TensorGenerator(...); select_undeuterated(T1.all_tp_clusters)
        n_undeut_runs = number of undeuterated HDX runs included in the 'library_info' master csv
        """

        if all_tp_clusters is None:
            all_tp_clusters = self.all_tp_clusters

        if name is None:
            name = self.name

        if library_info is None:
            library_info = self.library_info

        if n_undeut_runs is None:
            n_undeut_runs = self.n_undeut_runs

        my_seq = library_info.loc[library_info["name"] == name]["sequence"].values[0]

        if (
                self.old_data_dir is not None
        ):  # if comparing to old data, save old-data's fits in-place TODO: CONSIDER OUTPUTTING TO SNAKEMAKE DIR
            # open first three (undeut) dicts in list, store fit to theoretical dist
            for charge_dict in self.old_data:
                undeut_amds = [{
                    "major_species_integrated_intensities":
                        charge_dict["major_species_integrated_intensities"][i]
                } for i in range(3)]  # hardcode for gabe's undeut idxs in list
                charge_dict["fit_to_theo_dist"] = max(
                    [self.calculate_isotope_dist_dot_product(sequence=my_seq, undeut_integrated_mz_array=d) for d in undeut_amds]
                    )

        undeuts = []
        for ic in all_tp_clusters[
                0]:  # anticipates all undeuterated replicates being in the 0th index
            undeuts.append(ic)
        dot_products = []

        # <ake list of all normed dot products between an undeuterated IC and the theoretical distribution.
        for ic in undeuts:
            df = pd.DataFrame(
                ic.baseline_integrated_mz,
                columns=["major_species_integrated_intensities"],
            )
            fit = self.calculate_isotope_dist_dot_product(sequence=my_seq, undeut_integrated_mz_array=ic.baseline_integrated_mz)
            ic.undeut_ground_dot_product = fit
            dot_products.append((fit, ic.charge_states))

        # Append final (0, 0) to be called by charge_idxs which are not in the charge group for a single loop iteration
        dot_products.append((0, 0))
        charges = list(set(np.concatenate([ic.charge_states for ic in undeuts
                                          ])))
        out = dict.fromkeys(charges)
        charge_fits = dict.fromkeys(charges)
        for charge in charges:
            # print(charge)
            # Create sublist of undeuts with single charge state and same shape as undeuts, use -1 for non-matches to retain shape of list
            # [charge] == dot_products[i][1] ensures we only pick undeut_grounds from unconcatenated DataTensors, this saves trouble in undeut comparisons
            charge_idxs = []
            for i in range(len(dot_products)):
                if [charge] == dot_products[i][1]:
                    charge_idxs.append(i)
                else:
                    charge_idxs.append(-1)
            # print(charge_idxs)
            # print(np.asarray(dot_products)[charge_idxs])

            # Select best fit of charge state, append to output
            best = undeuts[charge_idxs[np.argmax(
                np.asarray(dot_products)[charge_idxs][:, 0])]]

            out[charge] = best
            charge_fits[charge] = max(
                np.asarray(dot_products)[charge_idxs][:, 0])

        self.undeut_grounds = out
        self.undeut_ground_dot_products = charge_fits


    def gaussian_function(self, x, H, A, x0, sigma):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


    def gauss_fit(self, x, y):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        mean = sum(x * y) / sum(y)
        sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
        nonzeros = [index for index, value in enumerate(list(y)) if value != 0]
        popt, pcov = curve_fit(self.gaussian_function, x, y, p0=[0, max(y), mean, sigma],
                               bounds=([0, 0, nonzeros[0], 0], [np.inf, np.inf, nonzeros[-1], np.inf]))
        return popt


    def rmse_from_gaussian_fit(self, distribution):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        try:
            xdata = [i for i in range(len(distribution))]
            ydata = distribution
            y_gaussian_fit = self.gaussian_function(xdata, *self.gauss_fit(xdata, ydata))
            rmse = mean_squared_error(ydata / max(ydata), y_gaussian_fit / max(y_gaussian_fit), squared=False)
            return rmse
        except:
            return 100


    def precalculate_fit_to_ground(self,
                                   all_tp_clusters=None,
                                   undeut_grounds=None):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        if all_tp_clusters is None:
            all_tp_clusters = self.all_tp_clusters
        if undeut_grounds is None:
            undeut_grounds = self.undeut_grounds

        for timepoint in all_tp_clusters:
            for ic in timepoint:
                undeut = undeut_grounds[ic.charge_states[0]]

                ic.dt_ground_err = abs(ic.dt_coms - undeut.dt_coms)
                ic.rt_ground_err = abs(ic.rt_com - undeut.rt_com)
                ic.auc_ground_err = ic.log_baseline_auc - undeut.log_baseline_auc
                ic.dt_ground_fit = max(
                    np.correlate(undeut.dt_norms[0], ic.dt_norms[0], mode='full'))
                ic.rt_ground_fit = max(np.correlate(undeut.rt_norm, ic.rt_norm, mode='full'))

                # these are pre calculated in the factor class
                # ic.dt_gaussian_rmse = self.rmse_from_gaussian_fit(ic.dt_norms[0])
                # ic.rt_gaussian_rmse = self.rmse_from_gaussian_fit(ic.rt_norm)

                ic.log_baseline_auc_diff = ic.log_baseline_auc - undeut.log_baseline_auc


    def generate_sample_paths(self):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        starts = np.linspace(0, 0.7, 8)
        sample_paths = []
        for start in starts:
            slopes = np.logspace(-0.5, 1.8, int(16 * (1.0 - start)))
            for slope in slopes:
                sample_paths.append(self.clusters_close_to_line(start, slope))

        self.sample_paths = [list(path) for path in set(sample_paths)]


    def clusters_close_to_line(
        self,
        start,
        slope,
        undeut_grounds=None,
        prefiltered_ics=None,
        max_peak_center=None,
    ):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        if undeut_grounds is None:
            undeut_grounds = self.undeut_grounds
        if prefiltered_ics is None:
            prefiltered_ics = self.prefiltered_ics
        if max_peak_center is None:
            max_peak_center = self.max_peak_center

        path = []
        uds = []
        fits = []
        # pick best-fitting undeut charge-state for all paths
        for key in undeut_grounds.keys():
            uds.append(undeut_grounds[key])
            fits.append(undeut_grounds[key].undeut_ground_dot_product)
        path.append(uds[np.argmax(fits)])

        # Relies on use of prefiltered_ics naming convention, when prefilter = 0, all_tp_clusters should have the first n_undeut_runs collapsed into a single list and be named and passed as prefiltered_ics
        xs = np.arange(len(prefiltered_ics))
        expected_centers = 3 + ((start + (1.0 - start) *
                                 (1 - np.exp(-xs * slope / len(xs)))) *
                                (max_peak_center - 3) * 0.85)

        # prefiltered_ics always has length n_deut_runs+1, all undeuterated are collapsed into PO.prefiltered_ics[0]
        for tp in range(1, len(xs)):
            try:
                peak_dist_rankings = []
                for ic in prefiltered_ics[tp]:
                    if ic.baseline_integrated_mz_com > max_peak_center:
                        continue
                    peak_dist_rankings.append((
                        (abs(ic.baseline_integrated_mz_com -
                             expected_centers[tp])),
                        ic,
                    ))
                peak_dist_rankings = sorted(peak_dist_rankings,
                                            key=lambda x: x[0])
                path.append(peak_dist_rankings[0][1])
            except:
                if len(peak_dist_rankings) == 0:
                    l = len(prefiltered_ics[tp])
                    if l > 0:
                        if l == 1:
                            # pick 0
                            path.append(prefiltered_ics[tp][0])
                        else:
                            # pick lowest mz
                            path.append(
                                sorted(
                                    prefiltered_ics[tp],
                                    key=lambda ic: ic.
                                    baseline_integrated_mz_com,
                                )[0])
                    else:
                        # No ics in tp, print message and append path[-1]
                        import os
                        import sys

                        # print("len(PO.prefiltered_ics["+str(tp)+"]) == 0")
                        path.append(path[-1])

        return tuple(path)


    def optimize_paths_multi(self, sample_paths=None, prefiltered_ics=None):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        # Main function of PO, returns the best-scoring HDX IC time-series 'path' of a set of bootstrapped paths.

        if sample_paths is None:
            sample_paths = self.sample_paths
        if prefiltered_ics is None:
            prefiltered_ics = self.prefiltered_ics

        final_paths = []
        for sample in sample_paths:
            current = copy.copy(sample)
            edited = True
            while edited:

                edited = False
                n_changes = 0
                ic_indices = []
                alt_paths = []

                for tp in range(1, len(current)):
                    for ic in prefiltered_ics[tp]:
                        buffr = copy.copy(current)
                        buffr[tp] = ic
                        alt_paths.append(buffr)

                # Decorate alt_paths
                combo_scoring = []
                for path in alt_paths:
                    combo_scoring.append(self.combo_score_multi(path))

                if min(combo_scoring) < self.combo_score_multi(current):
                    current = alt_paths[combo_scoring.index(min(combo_scoring))]
                    n_changes += 1
                    edited = True

                current_score = self.combo_score_multi(current)

                if edited == False:
                    final_paths.append(current)

        final_scores = []
        for path in final_paths:
            final_scores.append(self.combo_score_multi(path))

        # This order must be maintained, self.winner must exist before calling find_runners; winner and runners are both needed for set_bokeh tuple
        self.winner = final_paths[final_scores.index(min(final_scores))]
        self.winner_scores = self.report_score_mutli(self.winner)
        self.find_runners_multi()
        self.set_bokeh_tuples()
        self.filter_runners()
        # Compute the coefficient of variation (std.dev./mean) for RT and DT dimensions.
        self.rt_com_cv = (np.var([ic.rt_com for ic in self.winner if ic.rt_com is not None]) **
                          0.5) / np.mean([ic.rt_com for ic in self.winner if ic.rt_com is not None])
        self.dt_com_cv = (np.var([
            np.mean(ic.dt_coms) for ic in self.winner if ic.dt_coms is not None
        ])**0.5) / np.mean([np.mean(ic.dt_coms) for ic in self.winner if ic.dt_coms is not None])

    def find_runners_multi(self):
        # sets self.runners atr. sorts 'runner-up' single substitutions for each tp by score, lower is better.
        winner = self.winner
        prefiltered_ics = self.prefiltered_ics

        runners = []
        for tp in range(len(winner)):

            alt_paths = []
            for ic in prefiltered_ics[tp]:
                if ic is not winner[tp]:
                    buffr = copy.copy(winner)
                    buffr[tp] = ic
                    alt_paths.append(buffr)

            combo_scoring = []
            for path in alt_paths:
                combo_scoring.append(self.combo_score_multi(path))

            out_buffer = []
            for i in range(len(combo_scoring)):
                min_idx = combo_scoring.index(min(combo_scoring))
                out_buffer.append(alt_paths[min_idx][tp])
                alt_paths.pop(min_idx)
                combo_scoring.pop(min_idx)
            runners.append(out_buffer)

        self.runners = runners

    def optimize_paths_mono(self, sample_paths=None, prefiltered_ics=None):
        # Main function of PO, returns the best-scoring HDX IC time-series 'path' of a set of bootstrapped paths.

        if sample_paths is None:
            sample_paths = self.sample_paths
        if prefiltered_ics is None:
            prefiltered_ics = self.prefiltered_ics

        final_paths = []
        for sample in sample_paths[:1]:
            current = copy.copy(sample)
            edited = True
            while edited:

                edited = False
                n_changes = 0
                ic_indices = []
                alt_paths = []

                for tp in range(1, len(current)):
                    for ic in prefiltered_ics[tp]:
                        buffr = copy.copy(current)
                        buffr[tp] = ic
                        alt_paths.append(buffr)

                # Decorate alt_paths
                combo_scoring = []
                for pth in alt_paths:
                    combo_scoring.append(self.combo_score_mono(pth))

                if min(combo_scoring) < self.combo_score_mono(current):
                    current = alt_paths[combo_scoring.index(min(combo_scoring))]
                    n_changes += 1
                    edited = True

                current_score = self.combo_score_mono(current)

                if edited == False:
                    final_paths.append(current)
        final_scores = []
        for pth in final_paths:
            final_scores.append(self.combo_score_mono(pth))

        # This order must be maintained, self.winner must exist before calling find_runners; winner and runners are both needed for set_bokeh tuple
        self.winner = final_paths[final_scores.index(min(final_scores))]
        self.winner_scores = self.report_score_mono(self.winner)
        self.find_runners_mono()
        self.set_bokeh_tuples()
        self.filter_runners()
        self.rt_com_cv = (np.var([ic.rt_com for ic in self.winner if ic.rt_com is not None]) **
                          0.5) / np.mean([ic.rt_com for ic in self.winner if ic.rt_com is not None])
        self.dt_com_cv = (np.var([
            np.mean(ic.dt_coms) for ic in self.winner if ic.dt_coms is not None
        ]) ** 0.5) / np.mean([np.mean(ic.dt_coms) for ic in self.winner if ic.dt_coms is not None])
        # Doesn't return, only sets PO attributes

    def find_runners_mono(self):
        # sets self.runners atr. sorts 'runner-up' single substitutions for each tp by score, lower is better.
        winner = self.winner
        prefiltered_ics = self.prefiltered_ics

        runners = []
        for tp in range(len(winner)):

            alt_paths = []
            for ic in prefiltered_ics[tp]:
                if ic is not winner[tp]:
                    buffr = copy.copy(winner)
                    buffr[tp] = ic
                    alt_paths.append(buffr)

            combo_scoring = []
            for pth in alt_paths:
                combo_scoring.append(self.combo_score_mono(pth))

            out_buffer = []
            for i in range(len(combo_scoring)):
                min_idx = combo_scoring.index(min(combo_scoring))
                out_buffer.append(alt_paths[min_idx][tp])
                alt_paths.pop(min_idx)
                combo_scoring.pop(min_idx)
            runners.append(out_buffer)

        self.runners = runners

    def set_bokeh_tuples(self):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        # Sets IC.bokeh_tuple to be passed to bokeh for display through the HoverTool
        # Winners store the full values of the winning series scores
        # Runners store the differences between the winning scores and the score if they were to be substituted

        def score_dict(series):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            return {
                "int_mz_std_rmse":
                    self.int_mz_std_rmse(series) * self.int_mz_std_rmse_weight,
                "delta_mz_rate_backward":
                    self.delta_mz_rate(series)[0] * self.delta_mz_rate_backward_weight,
                "delta_mz_rate_afterward":
                    self.delta_mz_rate(series)[1] * self.delta_mz_rate_forward_weight,
                "dt_ground_rmse":
                    self.dt_ground_rmse(series) * self.dt_ground_rmse_weight,
                "rt_ground_rmse":
                    self.rt_ground_rmse(series) * self.rt_ground_rmse_weight,
                "dt_ground_fit":
                    self.dt_ground_fit(series) * self.dt_ground_fit_weight,
                "rt_ground_fit":
                    self.rt_ground_fit(series) * self.rt_ground_fit_weight,
                "baseline_peak_error":
                    self.baseline_peak_error(series) *
                    self.baseline_peak_error_weight,
                "auc_ground_rmse":
                    self.auc_ground_rmse(series) * self.auc_ground_rmse_weight,
                "rmses_sum":
                    self.rmses_sum(series) * self.rmses_sum_weight,
                "maxint_sum":
                    self.maxint_sum(series) * self.rmses_sum_weight,
                "int_mz_FWHM_rmse":
                    self.int_mz_FWHM_rmse(series) * self.int_mz_FWHM_rmse_weight,
                "nearest_neighbor_penalty":
                    self.nearest_neighbor_penalty(series) * self.nearest_neighbor_penalty_weight,
            }

        def score_diff(winner_scores, substituted_scores):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            return (
                winner_scores["int_mz_std_rmse"] -
                substituted_scores["int_mz_std_rmse"],
                winner_scores["delta_mz_rate_backward"] -
                substituted_scores["delta_mz_rate_backward"],
                winner_scores["delta_mz_rate_afterward"] -
                substituted_scores["delta_mz_rate_afterward"],
                winner_scores["dt_ground_rmse"] -
                substituted_scores["dt_ground_rmse"],
                winner_scores["rt_ground_rmse"] -
                substituted_scores["rt_ground_rmse"],
                winner_scores["dt_ground_fit"] -
                substituted_scores["dt_ground_fit"],
                winner_scores["rt_ground_fit"] -
                substituted_scores["rt_ground_fit"],
                winner_scores["baseline_peak_error"] -
                substituted_scores["baseline_peak_error"],
                winner_scores["auc_ground_rmse"] -
                substituted_scores["auc_ground_rmse"],
                winner_scores["rmses_sum"]
                - substituted_scores["rmses_sum"],
                winner_scores["maxint_sum"]
                - substituted_scores["maxint_sum"],
                winner_scores["int_mz_FWHM_rmse"]
                - substituted_scores["int_mz_FWHM_rmse"],
                winner_scores["nearest_neighbor_penalty"]
                - substituted_scores["nearest_neighbor_penalty"],



                sum([winner_scores[key] for key in winner_scores.keys()]) -
                sum([
                    substituted_scores[key]
                    for key in substituted_scores.keys()
                ]),
            )

        winner = self.winner
        runners = self.runners
        winner_scores = score_dict(winner)

        # Winners store absolute values of scores
        for ic in winner:
            ic.bokeh_tuple = ic.info_tuple + (
                ic.rt_ground_err,
                ic.dt_ground_err,
                winner_scores["int_mz_std_rmse"],
                winner_scores["delta_mz_rate_backward"],
                winner_scores["delta_mz_rate_afterward"],
                winner_scores["dt_ground_rmse"],
                winner_scores["rt_ground_rmse"],
                winner_scores["dt_ground_fit"],
                winner_scores["rt_ground_fit"],
                winner_scores["baseline_peak_error"],
                winner_scores["auc_ground_rmse"],
                winner_scores["rmses_sum"],
                winner_scores["maxint_sum"],
                winner_scores["int_mz_FWHM_rmse"],
                winner_scores["nearest_neighbor_penalty"],
                0,
            )

        # Runners store the differences between the winning scores and their scores when substituted
        # Negative values here mean the Runner causes that score to be worse, as all scores are low-is-better and the diffs are calculated winner-substituted
        for tp in range(len(runners)):
            for ic in runners[tp]:
                substituted_series = copy.copy(winner)
                substituted_series[tp] = ic
                substituted_scores = score_dict(substituted_series)
                ic.bokeh_tuple = (ic.info_tuple +
                                  (ic.rt_ground_err, ic.dt_ground_err) +
                                  score_diff(winner_scores, substituted_scores))

    def filter_runners(self, n_runners=5):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        filtered_runners = []
        for tp in self.runners:
            if len(tp) > n_runners:
                # Sum score diffs and sort index array in descending order by score value
                # Put first n_runners ics into filtered_runners
                tp_scores = [sum(ic.bokeh_tuple[-8:]) for ic in tp]
                tp_idxs = list(range(len(tp)))
                hi_to_lo = sorted(tp_idxs,
                                  key=lambda idx: tp_scores[idx],
                                  reverse=True)
                filtered_runners.append([
                    tp[idx]
                    for idx in hi_to_lo
                    if hi_to_lo.index(idx) < n_runners
                ])
            else:
                # if len(tp)<n_runners just append tp to filtered_runners
                filtered_runners.append(tp)
        self.filtered_runners = filtered_runners

    def calculate_theoretical_isotope_dist_from_sequence(self, sequence, n_isotopes=None):
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

    def calculate_empirical_isotope_dist_from_integrated_mz(self, integrated_mz_array,
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

    def calculate_isotope_dist_dot_product(self, sequence, undeut_integrated_mz_array):
        """Calculate dot product between theoretical isotope distribution from the sequence and experimental integrated mz array.
        
        Args:
            sequence (string): single-letter sequence of the library protein-of-interest
            undeut_integrated_mz_array (Numpy ndarray): observed integrated mz array from an undeuterated .mzML
        Returns:
            dot_product (float): result of dot product between theoretical and observed integrated-m/Z, from [0-1]

        """
        theo_isotope_dist = self.calculate_theoretical_isotope_dist_from_sequence(
            sequence=sequence)
        emp_isotope_dist = self.calculate_empirical_isotope_dist_from_integrated_mz(
            integrated_mz_array=undeut_integrated_mz_array)
        min_length = min([len(theo_isotope_dist), len(emp_isotope_dist)])
        dot_product = np.linalg.norm(
            np.dot(theo_isotope_dist[0:min_length], emp_isotope_dist[0:min_length])
        ) / np.linalg.norm(theo_isotope_dist) / np.linalg.norm(emp_isotope_dist)
        return dot_product

    ##########################################################################################################################################################################################################################################
    ### Scoring Functions for PathOptimizer ##################################################################################################################################################################################################
    ##########################################################################################################################################################################################################################################

    def int_mz_std_rmse(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        # calculates the difference in standard deviation from the mean from timepoint i-1 to i for i in [2, len(ics)]
        sd = 0
        for i in range(2, len(ics)):
            sd += (ics[i].baseline_integrated_mz_std -
                   ics[i - 1].baseline_integrated_mz_std)**2.0

        return math.sqrt(sd)

    def gabe_delta_mz_rate(self, major_species_centroids, timepoints=None):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        # reproduce logic of delta_mz_rate for Gabe's old data

        if timepoints is None:
            timepoints = self.timepoints

        # take mean of undeuts, fix for gabe's data
        major_species_centroids[0] = np.mean(
            [major_species_centroids.pop(0) for i in range(3)])
        sd = 0
        previous_rate = (major_species_centroids[1] - major_species_centroids[0]
                        ) / (timepoints[1] - timepoints[0])
        for i in range(2, len(major_species_centroids)):
            # if previous_rate == 0: diagnostic for /0 error
            new_com = major_species_centroids[i]
            if new_com < major_species_centroids[i - 1]:  # if we went backwards
                sd += (100 * (new_com - major_species_centroids[i - 1])**2.0
                      )  # penalize for going backwards
                new_com = (
                    major_species_centroids[i - 1] + 0.01
                )  # pretend we went forwards for calculating current rate
            current_rate = max([(new_com - major_species_centroids[i - 1]), 0.1
                               ]) / (timepoints[i] - timepoints[i - 1])
            if (current_rate / previous_rate) > 1.2:
                sd += (current_rate / previous_rate)**2.0
            previous_rate = current_rate
        return sd / len(major_species_centroids)

    def delta_mz_rate(self, ics, timepoints=None):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
    # Two penalizations are computed: [0] if the ic is too fast (sd) and [1] if the ic goes backwards (back)

        if timepoints is None:
            timepoints = self.timepoints
        
        backward = 0
        forward = 0
        previous_rate = max([(ics[1].baseline_integrated_mz_com - ics[0].baseline_integrated_mz_com) / (timepoints[1] - timepoints[0]), 0.1])

        for i in range(2, len(ics)):
            # if previous_rate == 0: diagnostic for /0 error
            new_com = ics[i].baseline_integrated_mz_com
            if new_com < ics[
                    i - 1].baseline_integrated_mz_com:  # if we went backwards
                backward += (100 *
                       (new_com - ics[i - 1].baseline_integrated_mz_com)**2.0
                      )  # penalize for going backwards
                new_com = (
                    ics[i - 1].baseline_integrated_mz_com + 0.01
                )  # pretend we went forwards for calculating current rate
            current_rate = max([
                (new_com - ics[i - 1].baseline_integrated_mz_com), 0.1
            ]) / (timepoints[i] - timepoints[i - 1])
            if (current_rate / previous_rate) > 1.2:
                forward += (current_rate / previous_rate)**2.0
            previous_rate = current_rate
        return backward / len(ics), forward / len(ics),

    def int_mz_rot_fit(self, ics):
        # Compares i to i-1 from ics[2]
        errors = []
        for i in range(2, len(ics)):
            i_mz, j_mz = (
                ics[i].baseline_integrated_mz_norm,
                ics[i - 1].baseline_integrated_mz_norm,
            )

            new_indices_i = np.nonzero(i_mz)[0] - np.argmax(i_mz)
            new_indices_j = np.nonzero(j_mz)[0] - np.argmax(j_mz)

            concat_indices = np.concatenate([new_indices_i, new_indices_j])
            common_low_index = min(concat_indices)
            common_high_index = max(concat_indices)

            new_array_i, new_array_j = (
                np.zeros((common_high_index - common_low_index + 1)),
                np.zeros((common_high_index - common_low_index + 1)),
            )

            new_indices_i -= common_low_index
            new_indices_j -= common_low_index

            new_array_i[new_indices_i] = i_mz[np.nonzero(i_mz)[0]]
            new_array_j[new_indices_j] = j_mz[np.nonzero(j_mz)[0]]

            errors.append(np.dot(new_array_i, new_array_j))
        return -np.average(errors)

    def dt_ground_rmse(
        self, ics
    ):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
      # rmse penalizes strong single outliers, score is minimized - lower is better
        return math.sqrt(sum([ic.dt_ground_err**2 for ic in ics]) / len(ics))

    def rt_ground_rmse(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        return math.sqrt(sum([ic.rt_ground_err**2 for ic in ics]) / len(ics))

    def dt_ground_fit(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        return sum([(1.0 / ic.dt_ground_fit) for ic in ics])

    def rt_ground_fit(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        return sum([(1.0 / ic.rt_ground_fit) for ic in ics])

    def baseline_peak_error(self, ics):  # Use RMSE instead TODO
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        # returns avg of peak_errors from baseline subtracted int_mz -> minimize score
        return np.average([ic.baseline_peak_error for ic in ics])

    def auc_ground_rmse(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        return np.sqrt(np.mean([ic.auc_ground_err**2 for ic in ics]))

    def auc_rmse(self,ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        sd = 0
        for ic in ics:
            sd += ic.log_baseline_auc_diff ** 2
        return math.sqrt(np.mean(sd))
    
    def rmses_sum(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        rmses = 0
        for ic in ics:
            rmses += 100*ic.baseline_integrated_mz_rmse
        return rmses
  
    def maxint_sum(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        maxint = 0
        for ic in ics:
            maxint += max(ic.baseline_integrated_mz)
            return 100000/maxint
  
    def int_mz_FWHM_rmse(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        sd = 0
        for i in range(2, len(ics)):
            sd += (
                ics[i].baseline_integrated_mz_FWHM - ics[i - 1].baseline_integrated_mz_FWHM
            ) ** 2.0
 
        return math.sqrt(sd)

    def nearest_neighbor_penalty(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        nn_penalty = 0
        for ic in ics:
            nn_penalty += 100 * (
                np.min([abs(1.0 - ic.nearest_neighbor_correlation), 0.5])
            ) ** 2.0
        return nn_penalty


    # Eventually put defaults here as else statements
    def set_score_weights(
        self,
        int_mz_std_rmse_weight=None,
        baseline_peak_error_weight=None,
        delta_mz_rate_backward_weight=None,
        delta_mz_rate_forward_weight=None,
        int_mz_rot_fit_weight=None,
        dt_ground_fit_weight=None,
        rt_ground_rmse_weight=None,
        dt_ground_rmse_weight=None,
        auc_ground_rmse_weight=None,
        rmses_sum_weight=None,
        maxint_sum_weight=None,
        int_mz_FWHM_rmse_weight=None,
        nearest_neighbor_penalty_weight=None
    ):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """

        if int_mz_std_rmse_weight != None:
            self.int_mz_std_rmse_weight = int_mz_std_rmse_weight
        if baseline_peak_error_weight != None:
            self.baseline_peak_error_weight = baseline_peak_error_weight
        if delta_mz_rate_backward_weight != None:
            self.delta_mz_rate_backward_weight = delta_mz_rate_backward_weight
        if delta_mz_rate_forward_weight != None:
            self.delta_mz_rate_forward_weight = delta_mz_rate_forward_weight
        if int_mz_rot_fit_weight != None:
            self.int_mz_rot_fit_weight = int_mz_rot_fit_weight
        if rt_ground_fit_weight != None:
            self.rt_ground_fit_weight = rt_ground_fit_weight
        if dt_ground_fit_weight != None:
            self.dt_ground_fit_weight = dt_ground_fit_weight
        if rt_ground_rmse_weight != None:
            self.rt_ground_rmse_weight = rt_ground_rmse_weight
        if dt_ground_rmse_weight != None:
            self.dt_ground_rmse_weight = dt_ground_rmse_weight
        if auc_ground_rmse_weight != None:
            self.auc_ground_rmse_weight = auc_ground_rmse_weight
        if rmses_sum_weight != None:
            self.rmses_sum_weight = rmses_sum_weight
        if maxint_sum_weight != None:
            self.maxint_sum_weight = maxint_sum_weight
        if int_mz_FWHM_rmse_weight != None:
            self.int_mz_FWHM_rmse_weight  = int_mz_FWHM_rmse_weight
        if nearest_neighbor_penalty_weight != None:
            self.nearest_neighbor_penalty_weight = nearest_neighbor_penalty_weight

    def combo_score_multi(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
        coeffs = [1.0, 0.10008302970275278, 0.09845279850378537, 0.772117631651534,
                  0.5310720751993441, 0.26070374293926435, 0.38586818508608206, 0.2522380132144197,
                  0.24236622286514023, 0.07221833422901867, 0.15051371151603132]

        return sum([
            coeffs[0] * self.baseline_peak_error_weight * self.baseline_peak_error(ics),
            coeffs[1] * self.delta_mz_rate_backward_weight * self.delta_mz_rate(ics)[0],
            coeffs[2] * self.delta_mz_rate_forward_weight * self.delta_mz_rate(ics)[1],
            coeffs[3] * self.dt_ground_rmse_weight * self.dt_ground_rmse(ics),
            coeffs[4] * self.dt_ground_fit_weight * self.dt_ground_fit(ics),
            coeffs[5] * self.rt_ground_fit_weight * self.rt_ground_fit(ics),
            coeffs[6] * self.rt_ground_rmse_weight * self.rt_ground_rmse(ics),
            coeffs[7] * self.auc_ground_rmse_weight * self.auc_ground_rmse(ics),
            coeffs[8] * self.rmses_sum_weight * self.rmses_sum(ics),
            coeffs[9] * self.int_mz_FWHM_rmse_weight * self.int_mz_FWHM_rmse(ics),
            coeffs[10] * self.nearest_neighbor_penalty_weight * self.nearest_neighbor_penalty(ics)
        ])

    def combo_score_mono(self, ics):
        coeffs = [1.0, 0.10008302970275278, 0.09845279850378537, 0.772117631651534,
                  0.5310720751993441, 0.26070374293926435, 0.38586818508608206, 0.2522380132144197,
                  0.24236622286514023, 0.07221833422901867, 0.15051371151603132]

        return sum([
            coeffs[0] * self.baseline_peak_error_weight * self.baseline_peak_error(ics),
            #coeffs[1] * self.delta_mz_rate_backward_weight * self.delta_mz_rate(ics)[0],
            #coeffs[2] * self.delta_mz_rate_forward_weight * self.delta_mz_rate(ics)[1],
            coeffs[3] * self.dt_ground_rmse_weight * self.dt_ground_rmse(ics),
            coeffs[4] * self.dt_ground_fit_weight * self.dt_ground_fit(ics),
            coeffs[5] * self.rt_ground_fit_weight * self.rt_ground_fit(ics),
            coeffs[6] * self.rt_ground_rmse_weight * self.rt_ground_rmse(ics),
            coeffs[7] * self.auc_ground_rmse_weight * self.auc_ground_rmse(ics),
            coeffs[8] * self.rmses_sum_weight * self.rmses_sum(ics),
            #coeffs[9] * self.int_mz_FWHM_rmse_weight * self.int_mz_FWHM_rmse(ics),
            coeffs[10] * self.nearest_neighbor_penalty_weight * self.nearest_neighbor_penalty(ics),
        ])
           
                                         
    def report_score_multi(self, ics):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """
    # TODO Add additional scores to this function                                    

        return {
            "int_mz_std_rmse": (self.int_mz_std_rmse_weight,
                                self.int_mz_std_rmse(ics)),
            "baseline_peak_error": (
                self.baseline_peak_error_weight,
                self.baseline_peak_error(ics),
            ),
            "delta_mz_rate_backard":
                (self.delta_mz_rate_backward_weight, self.delta_mz_rate(ics)[0]),
            "delta_mz_rate_foward":
                (self.delta_mz_rate_forward_weight, self.delta_mz_rate(ics)[1]),
            # self.int_mz_rot_fit_weight*self.int_mz_rot_fit(ics),
            "dt_ground_rmse": (self.dt_ground_rmse_weight,
                               self.dt_ground_rmse(ics)),
            "dt_ground_fit":
                (self.dt_ground_fit_weight, self.dt_ground_fit(ics)),
            "rt_ground_fit":
                (self.rt_ground_fit_weight, self.rt_ground_fit(ics)),
            "rt_ground_rmse": (self.rt_ground_rmse_weight,
                               self.rt_ground_rmse(ics)),
            "auc_ground_rmse": (self.auc_ground_rmse_weight,
                                self.auc_ground_rmse(ics)),
            "rmses_sum": (self.rmses_sum_weight,
                                self.rmses_sum(ics)),
            "int_mz_FWHM_rmse": (self.int_mz_FWHM_rmse_weight,
                                 self.int_mz_FWHM_rmse(ics)),
            "nearest_neighbor_penalty": (self.nearest_neighbor_penalty_weight,
                                         self.nearest_neighbor_penalty(ics)),
        }

    def report_score_mono(self, ics):
        # TODO Add additional scores to this function

        return {
            "baseline_peak_error": (
                self.baseline_peak_error_weight,
                self.baseline_peak_error(ics),
            ),
            "dt_ground_rmse": (self.dt_ground_rmse_weight,
                               self.dt_ground_rmse(ics)),
            "dt_ground_fit":
                (self.dt_ground_fit_weight, self.dt_ground_fit(ics)),
            "rt_ground_fit":
                (self.rt_ground_fit_weight, self.rt_ground_fit(ics)),
            "rt_ground_rmse": (self.rt_ground_rmse_weight,
                               self.rt_ground_rmse(ics)),
            "auc_ground_rmse": (self.auc_ground_rmse_weight,
                                self.auc_ground_rmse(ics)),
            "rmses_sum": (self.rmses_sum_weight,
                          self.rmses_sum(ics)),
            "int_mz_FWHM_rmse": (self.int_mz_FWHM_rmse_weight,
                          self.int_mz_FWHM_rmse(ics)),
            "nearest_neighbor_penalty": (self.nearest_neighbor_penalty_weight,
                                 self.nearest_neighbor_penalty(ics)),
        }

    def bokeh_plot(self, outpath):
        """Description of function.

        Args:
            arg_name (type): Description of input variable.

        Returns:
            out_name (type): Description of any returned objects.

        """

        def manual_cmap(value, low, high, palette):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            interval = (high - low) / len(palette)
            n_colors = len(palette)
            if value <= interval:
                return palette[0]
            else:
                if value > (n_colors - 1) * interval:
                    return palette[n_colors - 1]
                else:
                    for i in range(1, n_colors - 2):
                        if value > interval * i and value <= interval * i + 1:
                            return palette[i]

        def winner_added_mass_plotter(source, tooltips, old_source=None):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            p = figure(
                title=
                "Winning Timeseries Mean Added-Mass, Colored by RTxDT Error in ms",
                plot_height=400,
                plot_width=1275,
                background_fill_color="whitesmoke",
                x_range=(-1,
                         max([int(tp) for tp in source.data["timepoint"]]) + 1),
                tooltips=tooltips)
            err_mapper = linear_cmap(field_name="rtxdt_err",
                                     palette=Spectral6,
                                     low=0,
                                     high=1)
            color_bar = ColorBar(color_mapper=err_mapper["transform"],
                                 width=10,
                                 location=(0, 0))

            # Get mean value from source and map value to Spectral6
            mean_rtxdt_err = source.data["rtxdt_err"][0]
            mean_color = manual_cmap(mean_rtxdt_err, 0, 2, Spectral6)
            p.multi_line(
                xs="whisker_x",
                ys="whisker_y",
                source=source,
                line_color="black",
                line_width=1.5,
            )
            p.line(
                x="timepoint",
                y="baseline_integrated_mz_com",
                line_color=mean_color,
                source=source,
                line_width=3,
            )
            p.circle(
                x="timepoint",
                y="baseline_integrated_mz_com",
                source=source,
                line_color=err_mapper,
                color=err_mapper,
                fill_alpha=1,
                size=12,
            )

            if old_source is not None:  # plot added-masses of all charges of protein
                old_hover = HoverTool(
                    tooltips=[
                        ("Charge", "@charge"),
                        ("Delta MZ Rate Score", "@delta_mz_rate"),
                        (
                            "Fit of Undeuterated Added-Mass Distribution to Theoretical Distribution",
                            "@fit_to_theo_dist",
                        ),
                    ],
                    names=["old"],
                )

                old_ics = Line(
                    x="timepoint",
                    y="added_mass_centroid",
                    line_color="wheat",
                    line_width=1.5,
                )
                old_renderer = p.add_glyph(old_source, old_ics, name="old")
                p.add_tools(old_hover)

            p.xaxis.axis_label = "Timepoint Index"
            p.yaxis.axis_label = "Mean Added-Mass Units"
            p.min_border_top = 100
            p.min_border_left = 100
            p.min_border_right = 100
            p.add_layout(color_bar, "right")

            return p

        def winner_rtdt_plotter(source, tooltips):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            # set top margin
            p = figure(
                title=
                "Winning Timeseries RT and DT Center-of-Mass Error to Undeuterated Isotopic Cluster",
                plot_height=300,
                plot_width=1275,
                x_range=(-20, 20),
                y_range=(-1, 1),
                background_fill_color="whitesmoke",
                tooltips=tooltips)
            p.x(
                x="rt_ground_err",
                y="dt_ground_err",
                source=source,
                fill_alpha=1,
                size=5,
                color="black",
            )
            p.xaxis.axis_label = "RT Error (ms)"
            p.yaxis.axis_label = "DT Error (ms)"
            p.min_border_left = 100
            p.min_border_right = 100
            glyph = Text(x="rt_ground_err", y="dt_ground_err", text="timepoint")
            p.add_glyph(source, glyph)
            return p

        def winner_plotter(source, i, tooltips, old_source=None):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            if i == max([int(tp) for tp in source.data["timepoint"]]):
                p = figure(title="Timepoint " + str(i) +
                           ": Winning Isotopic-Cluster Added-Mass Distribution",
                           plot_height=400,
                           plot_width=450,
                           y_range=(0, 1),
                           background_fill_color="whitesmoke")
                p.min_border_bottom = 100
            else:
                p = figure(title="Timepoint " + str(i) +
                           ": Winning Isotopic Cluster Added-Mass Distribution",
                           plot_height=300,
                           plot_width=450,
                           y_range=(0, 1),
                           background_fill_color="whitesmoke")
            p.title.text_font_size = "8pt"
            index_view = CDSView(source=source,
                                 filters=[IndexFilter(indices=[i])])
            p.multi_line(xs="int_mz_x",
                         ys="int_mz_rescale",
                         source=source,
                         view=index_view,
                         line_color="blue",
                         line_width=1.5,
                         hover_color="red")
            p.add_tools(
                HoverTool(show_arrow=False,
                          line_policy="next",
                          tooltips=tooltips))

            # Have a figure by here, use glyph plotting from here
            """
            new_hover = HoverTool(tooltips=tooltips, names=["new"])
            index_view = CDSView(source=source, filters=[IndexFilter(indices=[i])])
            new_ics = MultiLine(
                xs="int_mz_x", ys="int_mz_rescale", line_color="blue", line_width=1.5
            )
            new_ics_hover = MultiLine(
                xs="int_mz_x", ys="int_mz_rescale", line_color="red", line_width=1.5
            )
            new_renderer = p.add_glyph(
                source, new_ics, view=index_view, name="new", hover_glyph=new_ics_hover
            )
            p.add_tools(new_hover)
            """

            if old_source is not None:  # plot ics matching the timepoint from old data
                old_hover = HoverTool(
                    tooltips=[
                        ("Charge", "@charge"),
                        ("Added-Mass Distribution Centroid", "@"),
                        ("Width", "@width"),
                    ],
                    names=["old"],
                )
                old_ics = MultiLine(
                    xs="int_mz_xs",
                    ys="int_mz_rescale",
                    line_color="wheat",
                    line_width=1.5,
                )
                old_ics_hover = MultiLine(
                    xs="int_mz_xs",
                    ys="int_mz_rescale",
                    line_color="red",
                    line_width=1.5,
                )
                old_tp_view = CDSView(
                    source=old_source,
                    filters=[
                        GroupFilter(column_name="type", group="ic"),
                        GroupFilter(column_name="timepoint", group=str(i)),
                    ],
                )
                old_renderer = p.add_glyph(
                    old_source,
                    old_ics,
                    view=old_tp_view,
                    hover_glyph=old_ics_hover,
                    name="old",
                )
                p.add_tools(old_hover)

            p.xaxis.axis_label = "Added-Mass Units"
            p.yaxis.axis_label = "Relative Intensity"
            p.min_border_left = 100
            return p

        def runner_plotter(source, i, tooltips):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            if i == max([int(tp) for tp in source.data["timepoint"]]):
                p = figure(
                    title="Runner-Up Isotopic Cluster Added-Mass Distributions",
                    plot_height=400,
                    plot_width=375,
                    y_range=(0, 1),
                    background_fill_color="whitesmoke")
                p.min_border_bottom = 100
            else:
                p = figure(
                    title="Runner-Up Isotopic Cluster Added-Mass Distributions",
                    plot_height=300,
                    plot_width=375,
                    y_range=(0, 1),
                    background_fill_color="whitesmoke",
                    tools="pan,wheel_zoom,hover,reset,help",
                    tooltips=tooltips)
            p.title.text_font_size = "8pt"
            runner_timepoint_view = CDSView(
                source=source,
                filters=[GroupFilter(column_name="timepoint", group=str(i))],
            )
            p.multi_line(
                xs="int_mz_x",
                ys="int_mz_rescale",
                source=source,
                view=runner_timepoint_view,
                line_color="blue",
                alpha=0.5,
                hover_color="red",
                hover_alpha=1,
                line_width=1.5,
            )
            p.add_tools(
                HoverTool(show_arrow=False,
                          line_policy="next",
                          tooltips=tooltips))
            p.xaxis.axis_label = "Added-Mass Units"
            p.yaxis.axis_label = "Relative Intensity"
            return p

        def rtdt_plotter(source, i, tooltips):
            """Description of function.

            Args:
                arg_name (type): Description of input variable.

            Returns:
                out_name (type): Description of any returned objects.

            """
            if i == max([int(tp) for tp in source.data["timepoint"]]):
                p = figure(title="RT and DT Error from Undeuterated",
                           plot_height=400,
                           plot_width=450,
                           background_fill_color="whitesmoke",
                           x_range=(-30, 30),
                           y_range=(-2, 2),
                           tooltips=tooltips)
                p.min_border_bottom = 100
            else:
                p = figure(
                    title=
                    "Retention and Drift Center-of-Mass Error to Undeuterated",
                    plot_height=300,
                    plot_width=450,
                    background_fill_color="whitesmoke",
                    x_range=(-30, 30),
                    y_range=(-2, 2),
                    tooltips=tooltips)
            p.title.text_font_size = "8pt"
            timepoint_runner_view = CDSView(
                source=source,
                filters=[
                    GroupFilter(column_name="timepoint", group=str(i)),
                    GroupFilter(column_name="winner_or_runner", group=str(1)),
                ],
            )
            p.circle(
                x="rt_ground_err",
                y="dt_ground_err",
                source=source,
                view=timepoint_runner_view,
                line_color="blue",
                hover_color="red",
                alpha=0.25,
                hover_alpha=1,
                size=5,
            )
            timepoint_winner_view = CDSView(
                source=source,
                filters=[
                    GroupFilter(column_name="timepoint", group=str(i)),
                    GroupFilter(column_name="winner_or_runner", group=str(0)),
                ],
            )
            p.circle(
                x="rt_ground_err",
                y="dt_ground_err",
                source=source,
                view=timepoint_winner_view,
                line_color="black",
                fill_color="black",
                hover_color="red",
                size=5,
            )
            p.xaxis.axis_label = "RT COM Error from Ground (ms)"
            p.yaxis.axis_label = "DT COM Error from Ground (ms)"
            p.min_border_right = 100
            return p

        output_file(outpath, mode="inline")

        # Start old_data source creation
        if self.old_data_dir is not None:

            # create bokeh datasource from gabe's hx_fits

            # divide source by data-types: single-tp ic-level data and charge-state time-series level data
            # old_charges will be used for plotting added-mass and time-series stats

            # init dicts with columns for plotting
            old_ics = dict.fromkeys([
                "timepoint",
                "added_mass_centroid",
                "added_mass_width",
                "int_mz_ys",
                "int_mz_xs",
                "type",
                "int_mz_rescale",
            ])
            for key in old_ics.keys():
                old_ics[key] = []

            old_charges = dict.fromkeys([
                "major_species_integrated_intensities",
                "centroid",
                "major_species_widths",
                "fit_to_theo_dist",
                "delta_mz_rate",
                "added_mass_xs",
                "lowers",
                "uppers",
                "type",
                "charge",
            ])
            for key in old_charges.keys():
                old_charges[key] = []

            # set switch to pull values that only need to be computed once, append each old_file's values to old_data{}
            old_switch = None
            for ts in self.old_data:

                int_mz_xs = list(
                    range(len(ts["major_species_integrated_intensities"][0])))
                timepoints = list(range(len(ts["major_species_centroid"])))
                print("old_data timepoints: " + str(timepoints))

                # Add line to old_charges for each charge file
                for key in [
                        "major_species_integrated_intensities",
                        "centroid",
                        "fit_to_theo_dist",
                        "delta_mz_rate",
                        "charge",
                ]:
                    old_charges[key].append(ts[key])

                old_charges["added_mass_xs"].append(timepoints)
                old_charges["major_species_widths"].append([
                    len(np.nonzero(ic)[0])
                    for ic in ts["major_species_integrated_intensities"]
                ])
                old_charges["lowers"].append([
                    ts["major_species_centroid"][tp] -
                    (ts["major_species_widths"][tp] / 2) for tp in timepoints
                ])
                old_charges["uppers"].append([
                    ts["major_species_centroid"][tp] +
                    (ts["major_species_widths"][tp] / 2) for tp in timepoints
                ])
                old_charges["type"].append("ts")

                # Add line to old_ics for each hdx timepoint in each charge file
                for tp in timepoints:
                    if tp < 3:
                        old_ics["timepoint"].append(str(0))
                    else:
                        old_ics["timepoint"].append(str(tp - 2))
                    old_ics["added_mass_centroid"].append(
                        ts["major_species_centroid"][tp])
                    old_ics["added_mass_width"].append(
                        len(
                            np.nonzero(
                                ts["major_species_integrated_intensities"][tp])
                            [0]))
                    old_ics["int_mz_ys"].append(
                        ts["major_species_integrated_intensities"][tp])
                    old_ics["int_mz_rescale"].append(
                        ts["major_species_integrated_intensities"][tp] /
                        max(ts["major_species_integrated_intensities"][tp]))
                    old_ics["int_mz_xs"].append(int_mz_xs)
                    old_ics["type"].append("ic")

            ts_df = pd.DataFrame.from_dict(
                old_charges
            )  # len = number of identified charge states for given protein name
            ic_df = pd.DataFrame.from_dict(
                old_ics)  # len = n_charges * n_hdx_timepoints

            self.old_undeut_ground_dot_products = dict.fromkeys(
                ts_df["charge"].values)
            for charge in self.old_undeut_ground_dot_products.keys():
                self.old_undeut_ground_dot_products[charge] = ts_df.loc[
                    ts_df["charge"] == charge]["fit_to_theo_dist"].values

            old_df = pd.concat([ts_df, ic_df])
            self.old_df = old_df
            # make cds from df
            gds = ColumnDataSource(old_df)

        else:
            old_df = pd.DataFrame()

        # End old_data source creation

        # TODO: Eventually all the info tuples should be dicts for ease of use in pandas and bokeh, change this style of source construction to dicts of lists
        # This has also just become a clusterfuck and needs to be cleaned up
        winner_data = []
        runner_data = []
        all_data = []

        winner_rtxdt_rmse = np.sqrt(
            np.mean([((ic.bokeh_tuple[18] * 0.07) * ic.bokeh_tuple[19])**2
                     for ic in self.winner]))
        for tp in range(len(self.winner)):
            edit_buffer = copy.copy(self.winner[tp].bokeh_tuple)
            edit_buffer = (edit_buffer[:18] + (edit_buffer[18] * 0.07,) +
                           edit_buffer[19:] + (
                               str(tp),
                               np.nonzero(edit_buffer[17])[0][-1],
                               np.nonzero(edit_buffer[17])[0][0],
                               "0",
                               ((edit_buffer[18] * 0.07) * edit_buffer[19]),
                               winner_rtxdt_rmse,
                               np.asarray([tp, tp]),
                               np.asarray([
                                   np.nonzero(edit_buffer[17])[0][0],
                                   np.nonzero(edit_buffer[17])[0][-1],
                               ]),
                               edit_buffer[17] / max(edit_buffer[17]),
                           ))  # 0.07 is adjustment from bins to ms
            winner_data.append(edit_buffer)
            all_data.append(edit_buffer)

        for tp in range(len(self.filtered_runners)):
            for ic in self.filtered_runners[tp]:
                edit_buffer = copy.copy(ic.bokeh_tuple)
                edit_buffer = (edit_buffer[:18] + (edit_buffer[18] * 0.07,) +
                               edit_buffer[19:] + (
                                   str(tp),
                                   np.nonzero(edit_buffer[17])[0][-1],
                                   np.nonzero(edit_buffer[17])[0][0],
                                   "1",
                                   ((edit_buffer[18] * 0.07) * edit_buffer[19]),
                                   "NA",
                                   np.asarray([tp, tp]),
                                   np.asarray([
                                       np.nonzero(edit_buffer[17])[0][0],
                                       np.nonzero(edit_buffer[17])[0][-1],
                                   ]),
                                   edit_buffer[17] / max(edit_buffer[17]),
                               ))  # 0.07 is adjustment from bins to ms
                runner_data.append(edit_buffer)
                all_data.append(edit_buffer)

        columns = [
            "source_file",
            "tensor_idx",
            "n_factors",
            "factor_idx",
            "cluster_idx",
            "charge_states",
            "n_concatenated",
            "mz_bin_low",  # from factor bins
            "mz_bin_high",  # from factor bins
            "baseline_subtracted_area_under_curve",
            "baseline_subtracted_grate_sum",
            "baseline_subtracted_peak_error",
            "baseline_integrated_mz_com",
            "abs_mz_com",
            "rt",
            "dt",
            "int_mz_x",
            "int_mz_y",
            "rt_ground_err",
            "dt_ground_err",
            "int_mz_std_rmse",
            "delta_mz_rate",
            "dt_ground_rmse_score",
            "rt_ground_rmse_score",
            "dt_ground_fit",
            "rt_ground_fit",
            "baseline_peak_error",
            "auc_ground_rmse",
            "net_score_difference",
            "timepoint",
            "upper_added_mass",
            "lower_added_mass",
            "winner_or_runner",
            "rtxdt_err",
            "rtxdt_rmse",
            "whisker_x",
            "whisker_y",
            "int_mz_rescale",
        ]

        winner_frame = pd.DataFrame(winner_data, columns=columns)
        runner_frame = pd.DataFrame(runner_data, columns=columns)
        all_frame = pd.DataFrame(all_data, columns=columns)

        wds = ColumnDataSource(winner_frame)
        rds = ColumnDataSource(runner_frame)
        ads = ColumnDataSource(all_frame)

        max_intensity = max(
            [max(int_mz) for int_mz in all_frame["int_mz_y"].values])

        # HoverToolTips, determines information to be displayed when hovering over a glyph
        winner_tts = [
            ("Tensor Index", "@tensor_idx"),
            ("Charge State(s)", "@charge_states"),
            ("Timepoint", "@timepoint"),
            ("Peak Error", "@baseline_subtracted_peak_error"),
            ("Center of Mass in Added-Mass_Units", "@baseline_integrated_mz_com"),
            ("Center of Mass in M/Z", "@abs_mz_com"),
            ("Retention Time COM Error to Ground", "@rt_ground_err"),
            ("Drift Time COM Error to Ground", "@dt_ground_err"),
            ("int_mz_std_rmse", "@int_mz_std_rmse"),
            ("delta_mz_rate", "@delta_mz_rate"),
            ("dt_ground_rmse", "@dt_ground_rmse_score"),
            ("rt_ground_rmse", "@rt_ground_rmse_score"),
            ("dt_ground_fit", "@dt_ground_fit"),
            ("rt_ground_fit", "@rt_ground_fit"),
            ("baseline_peak_error", "@baseline_peak_error"),
            ("auc_ground_rmse", "@auc_ground_rmse"),
            ("Δ scores", "positive is better"),
            ("Δ_net_score", "@net_score_difference"),
        ]

        runner_tts = [
            ("Tensor Index", "@tensor_idx"),
            ("Charge State(s)", "@charge_states"),
            ("Timepoint", "@timepoint"),
            ("Peak Error", "@baseline_subtracted_peak_error"),
            ("Center of Mass in Added-Mass_Units", "@baseline_integrated_mz_com"),
            ("Center of Mass in M/Z", "@abs_mz_com"),
            ("Retention Time COM Error to Ground", "@rt_ground_err"),
            ("Drift Time COM Error to Ground", "@dt_ground_err"),
            ("Δ scores", "positive is better"),
            ("Δ int_mz_std_err", "@int_mz_std_err"),
            ("Δ delta_mz_rate", "@delta_mz_rate"),
            ("Δ dt_ground_rmse", "@dt_ground_rmse_score"),
            ("Δ rt_ground_rmse", "@rt_ground_rmse_score"),
            ("Δ dt_ground_fit", "@dt_ground_fit"),
            ("Δ rt_ground_fit", "@rt_ground_fit"),
            ("Δ baseline_peak_error", "@baseline_peak_error"),
            ("Δ auc_ground_rmse", "@auc_ground_rmse"),
            ("Δ_net_score", "@net_score_difference"),
        ]

        mass_added_tts = [
            ("Timepoint", "@timepoint"),
            ("Charge State(s)", "@charge_states"),
            ("RT COM Error (ms)", "@rt_ground_err"),
            ("DT COM Error (ms)", "@dt_ground_err"),
            ("DTxRT Error", "@rtxdt_err"),
        ]

        winner_rtdt_tts = [
            ("Timepoint", "@timepoint"),
            ("Charge State(s)", "@charge_states"),
            ("RT COM Error (ms)", "@rt_ground_err"),
            ("DT COM Error (ms)", "@dt_ground_err"),
        ]

        n_timepoints = len(self.winner)
        # print("internal n_timepoints: "+str(n_timepoints))
        if self.old_data_dir is not None:
            winner_plots = [
                winner_plotter(wds, i, winner_tts, old_source=gds)
                for i in range(n_timepoints)
            ]

        else:
            winner_plots = [
                winner_plotter(wds, i, winner_tts) for i in range(n_timepoints)
            ]

        runner_plots = [
            runner_plotter(rds, i, runner_tts) for i in range(n_timepoints)
        ]
        rtdt_plots = [
            rtdt_plotter(ads, i, runner_tts) for i in range(n_timepoints)
        ]

        rows = []
        if self.old_data_dir is not None:
            rows.append(
                gridplot(
                    [
                        winner_added_mass_plotter(
                            wds, mass_added_tts, old_source=gds)
                    ],
                    sizing_mode="fixed",
                    toolbar_location="left",
                    ncols=1,
                ))

        else:
            rows.append(
                gridplot(
                    [winner_added_mass_plotter(wds, mass_added_tts)],
                    sizing_mode="fixed",
                    toolbar_location="left",
                    ncols=1,
                ))

        rows.append(
            gridplot(
                [winner_rtdt_plotter(wds, winner_rtdt_tts)],
                sizing_mode="fixed",
                toolbar_location="left",
                ncols=1,
            ))

        [
            rows.append(
                gridplot(
                    [winner_plots[i], runner_plots[i], rtdt_plots[i]],
                    sizing_mode="fixed",
                    toolbar_location="left",
                    ncols=3,
                )) for i in range(n_timepoints)
        ]

        if self.old_files is not None:
            final = column(
                Div(text=
                    """<h1 style='margin-left: 300px'>HDX Timeseries Plot for """
                    + self.name + """</h1>"""),
                Div(text=
                    "<h3 style='margin-left: 300px'>New Undeuterated-Ground Fits to Theoretical MZ Distribution: </h3>"
                   ),
                Div(text="<h3 style='margin-left: 300px'>" +
                    str(self.undeut_ground_dot_products) + "</h3>"),
                Div(text=
                    "<h3 style='margin-left: 300px'>Old Undeuterated-Ground Fits to Theoretical MZ Distribution: </h3>"
                   ),
                Div(text="<h3 style='margin-left: 300px'>" +
                    str(self.old_undeut_ground_dot_products) + "</h3>"),
                gridplot(rows,
                         sizing_mode="fixed",
                         toolbar_location=None,
                         ncols=1),
            )

        else:
            final = column(
                Div(text=
                    """<h1 style='margin-left: 300px'>HDX Timeseries Plot for """
                    + self.name + """</h1>"""),
                Div(text=
                    "<h3 style='margin-left: 300px'>Undeuterated-Ground Fits to Theoretical MZ Distribution: </h3>"
                   ),
                Div(text="<h3 style='margin-left: 300px'>" +
                    str(self.undeut_ground_dot_products) + "</h3>"),
                gridplot(rows,
                         sizing_mode="fixed",
                         toolbar_location=None,
                         ncols=1),
            )

        save(final)
