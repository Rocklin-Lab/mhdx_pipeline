import os
import time
import math
import copy
import psutil
import peakutils
import numpy as np
from nn_fac import ntf
from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.measurements import center_of_mass


class DataTensor:

    def __init__(self, source_file, tensor_idx, timepoint_idx, name,
                 total_mass_window, n_concatenated, charge_states, **kwargs):

        ###Set Common Attributes###

        # Positional args
        self.source_file = source_file
        self.tensor_idx = tensor_idx
        self.timepoint_idx = timepoint_idx
        self.name = name
        self.total_mass_window = total_mass_window
        self.n_concatenated = n_concatenated
        self.charge_states = charge_states
        self.measured_precision = 1e-6  # hardcode from pymzml
        self.internal_precision = (
            1000  # chosen for low-loss compression of mz dimension
        )

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
            if "lows" in kws:
                self.lows = kwargs["lows"]
            if "highs" in kws:
                self.highs = kwargs["highs"]
            if "abs_mz_low" in kws:
                self.mz_bin_low = kwargs["abs_mz_low"]

        ###Compute Instance Values###

        # Handle normal case: new DataTensor from output of isolate_tensors.py
        if self.n_concatenated == 1:
            """
            t0 = time.time()
            #print('Reprofiling...')
            self.grid_out = np.reshape(self.reprofile(), (len(self.rts), len(self.dts)))
            t = time.time()
            #print("T+"+str(t-t0))
            #For creating full_grid_out, a 3d array with all mz dimensions having the same length
            #Find min and max mz range values:
            #check high
            self.mz_bin_low = 1e9
            self.mz_bin_high = 0
            for x in self.grid_out:
                for y in x:
                    if len(y) > 0:
                        z = np.asarray(y)[:,0]
                        if min(z) < self.mz_bin_low: self.mz_bin_low = y[np.argmin(z)][0]
                        if max(z) > self.mz_bin_high: self.mz_bin_high = y[np.argmax(z)][0]

            #create zero array with range of bin indices
            self.mz_bins = np.arange(self.mz_bin_low, self.mz_bin_high, 0.02) #hardcode for desired coarseness of gaussian representation, could be parameterized
            self.mz_len = len(self.mz_bins)

            #create empty space with dimensions matching grid_out and m/z indices
            self.full_grid_out = np.zeros((np.shape(self.grid_out)[0], np.shape(self.grid_out)[1], self.mz_len))

            #determine and apply mapping of nonzero grid_out values to full_grid_out
            for l in range(len(self.grid_out)):
                for m in range(len(self.grid_out[l])):
                    if len(self.grid_out[l][m]) == 0:
                        pass

                    else:
                        if len(self.grid_out[l][m]) != 0:
                            low_idx = np.searchsorted(self.mz_bins, self.grid_out[l][m][0,0])
                            high_idx = np.searchsorted(self.mz_bins, self.grid_out[l][m][-1,0])
                            if high_idx - low_idx == len(self.grid_out[l][m][:,0]):
                                self.indices = np.arange(low_idx, high_idx, 1)
                            else:
                                self.indices = np.clip(np.searchsorted(self.mz_bins, self.grid_out[l][m][:,0]), 0, len(self.mz_bins)-1)
                        else:
                            self.indices = np.clip(np.searchsorted(self.mz_bins, self.grid_out[l][m][:,0]), 0, len(self.mz_bins)-1)

                        if self.indices[-1] == len(self.mz_bins):
                            self.indices[-1] = self.indices[-1]-1
                        try:
                            self.full_grid_out[l,m][self.indices] = self.grid_out[l][m][:,1]
                        except:
                            ipdb.set_trace()
                            print("this stops the iterator")
            """

            (
                self.retention_labels,
                self.drift_labels,
                self.mz_labels,
                self.min_mz,
                self.max_mz,
                self.full_grid_out,
            ) = self.sparse_to_full_tensor((self.rts, self.dts, self.seq_out))
            self.full_gauss_grids = self.gauss(self.full_grid_out)

        # Handle concatenated tensor case, check for required inputs
        else:
            if not all("dts" in kwargs and "rts" in kwargs and
                       "lows" in kwargs and "highs" in kwargs and
                       "concatenated_grid" in kwargs and
                       "abs_mz_low" in kwargs and "concat_dt_idxs" in kwargs):

                print("Concatenated Tensor Missing Required Values")
                sys.exit()

    def sparse_to_full_tensor(self, data):
        retention_labels, drift_labels, sparse_data = data

        min_mz = min([min(x[:, 0]) for x in sparse_data if len(x[:, 0]) > 0])
        max_mz = max([max(x[:, 0]) for x in sparse_data if len(x[:, 0]) > 0])
        # print (min_mz, max_mz)
        # mz_labels = np.arange(min_mz, max_mz + 0.03, 0.02)
        mz_labels = np.arange(min_mz, max_mz + 0.03, 0.02)

        tensor3_out = np.zeros(
            (len(retention_labels), len(drift_labels), len(mz_labels)))

        scan = 0
        for i in range(len(retention_labels)):
            for j in range(len(drift_labels)):
                mz_indices = np.searchsorted(mz_labels, sparse_data[scan][:, 0])
                large_array = np.zeros((len(sparse_data[scan]), len(mz_labels)))
                large_array[np.arange(len(sparse_data[scan])),
                            mz_indices] = sparse_data[scan][:, 1]
                tensor3_out[i][j] = np.sum(large_array, axis=0)
                scan += 1

        return (retention_labels, drift_labels, mz_labels, min_mz, max_mz,
                tensor3_out)

    def reprofile(self):
        # Reads list of mz peak centroids and returns full length mz array of gaussians.
        # Measured and interal precision are adjusted to vary compression of the mz dimension.
        # Adapted from pymzml, used here to push peak reprofiling downstream of mzml extraction for more efficient load distribution.
        out = []
        for scan in self.seq_out:
            tmp = ddict(int)
            for mz, i in scan:
                # Let the measured precision be 2 sigma of the signal width
                # When using normal distribution
                # FWHM = 2 sqt(2 * ln(2)) sigma = 2.3548 sigma
                s = mz * self.measured_precision * 2  # in before 2
                s2 = s * s
                floor = mz - 5.0 * s  # Gauss curve +- 3 sigma
                ceil = mz + 5.0 * s
                ip = self.internal_precision / 4
                # more spacing, i.e. less points describing the gauss curve
                # -> faster adding
                for _ in range(int(round(floor * ip)),
                               int(round(ceil * ip)) + 1):
                    if _ % int(5) == 0:
                        a = float(_) / float(ip)
                        y = i * math.exp(-1 * ((mz - a) * (mz - a)) / (2 * s2))
                        tmp[a] += y
            out.append(np.asarray([[key, tmp[key]] for key in list(tmp.keys())
                                  ]))
        return np.asarray(out)

    # Takes tensor input and gaussian filter parameters, outputs filtered data
    def gauss(self, grid, rt_sig=3, dt_sig=1):

        gauss_grid = np.zeros(np.shape(grid))
        for i in range(np.shape(grid)[2]):
            gauss_grid[:, :, i] = gaussian_filter(grid[:, :, i],
                                                  (rt_sig, dt_sig))
        return gauss_grid

    def interpolate(self, grid_in, new_mz_len, gauss_params=None):
        # Takes length of mz_bins to interpolated to, and optional gaussian filter parameters
        # Returns the interpolated tensor, length of interpolated axis, interpolated low_lims and high_lims

        if gauss_params != None:
            grid = self.gauss(grid_in, gauss_params[0], gauss_params[1])
        else:
            grid = grid_in

        test_points = []
        z_axis = np.clip(np.linspace(0,
                                     np.shape(grid)[2], new_mz_len), 0,
                         np.shape(grid)[2] - 1)
        for n in range(np.shape(grid)[0]):
            for o in range(np.shape(grid)[1]):
                for p in z_axis:
                    test_points.append((n, o, p))

        x, y, z = (
            np.arange(np.shape(grid)[0]),
            np.arange(np.shape(grid)[1]),
            np.arange(np.shape(grid)[2]),
        )
        interpolation_function = sp.interpolate.RegularGridInterpolator(
            points=[x, y, z], values=grid)
        interpolated_out = interpolation_function(test_points)
        interpolated_out = np.reshape(
            interpolated_out,
            (np.shape(grid)[0], np.shape(grid)[1], new_mz_len))

        interpolated_bin_mzs = np.linspace(self.mz_bin_low, self.mz_bin_high,
                                           new_mz_len)
        interpolated_low_lims = np.searchsorted(interpolated_bin_mzs,
                                                self.mz_labels[self.lows])
        interpolated_high_lims = np.searchsorted(interpolated_bin_mzs,
                                                 self.mz_labels[self.highs])

        return [interpolated_out, interpolated_low_lims, interpolated_high_lims]

    def factorize(self, n_factors=13, new_mz_len=None, gauss_params=None):
        # Test factorization starting at n_factors = 15 and counting down, keep factorization that has no factors with correlation greater than 0.2 in any dimension.

        def corr_check(factors, cutoff):
            # Checks scipy non_negatve_parafac output factors for inter-factor (off-diagonal) correlations > cutoff, returns True if all values are < cutoff

            a = np.minimum(
                np.minimum(np.corrcoef(factors[0].T),
                           np.corrcoef(factors[1].T)),
                np.corrcoef(factors[2].T),
            )

            if any(a[np.where(~np.eye(a.shape[0], dtype=bool))] > cutoff):
                return False
            else:
                return True

        def pmem(id_str):
            process = psutil.Process(os.getpid())
            print(id_str + " Process Memory (GB): " +
                  str(process.memory_info().rss / 1024 / 1024 / 1024))

        t = time.time()
        pmem("0 Start")
        # print('Filtering... T+'+str(t-t0))
        # handle concatenation and intetrpolfilter option
        if self.n_concatenated != 1:
            grid, lows, highs, concat_dt_idxs = (
                self.concatenated_grid,
                self.lows,
                self.highs,
                self.concat_dt_idxs,
            )
        else:
            if new_mz_len != None:
                if gauss_params != None:
                    grid, lows, highs, concat_dt_idxs = (
                        interpolate(
                            self.full_grid_out,
                            new_mz_len,
                            gauss_params[0],
                            gauss_params[1],
                        ),
                        None,
                    )
                else:
                    grid, lows, highs, concat_dt_idxs = (
                        interpolate(self.full_grid_out, new_mz_len),
                        None,
                    )
            else:
                lows, highs, concat_dt_idxs = self.lows, self.highs, None
                if gauss_params != None:
                    grid = self.gauss(self.full_grid_out, gauss_params[0],
                                      gauss_params[1])
                else:
                    grid = self.full_grid_out
        grid = self.full_gauss_grids
        pmem("1 Read Params")
        t = time.time()
        # print('Zeroing Non-POI M/z... T+'+str(t-t0))
        # Multiply all values outside of integration box boundaries by 0, TODO: demonstrate this against keeping the full tensor - obv faster, self-evidently better fac: quantify as support
        zero_mult = np.zeros((np.shape(grid)))
        for lo, hi in zip(lows, highs):
            zero_mult[:, :, lo:hi] = 1
        grid *= zero_mult
        pmem("2 Zeroing")
        # Count down from 15 and keep highest n_factors that satisfies corr_check
        flag = True
        t = time.time()
        # print('Start Factorization Series... T+'+str(t-t0))
        pmem("3 Pre-Factorization")
        n_itr = 4
        while flag:
            pmem(str(n_itr) + " " + str(n_factors) + " Factors " + " Start")
            t1 = time.time()
            # print('Starting '+str(nf)+' Factors... T+'+str(t1-t))
            nnf1 = ntf.ntf(grid, n_factors)
            pmem(str(n_itr) + " " + str(n_factors) + " Factors " + " End")
            n_itr += 1
            t2 = time.time()
            # print('Factorization Duration: '+str(t2-t1))

            if n_factors > 1:
                if corr_check(nnf1, 0.25):
                    flag = False
                    break
                else:
                    n_factors -= 1
            else:
                flag = False
                print("All n-factors failed for Index: " + str(self.name) +
                      ", keeping 1 factor decomposition.")
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
                    factor_idx=i,
                    n_factors=n_factors,
                    lows=lows,
                    highs=highs,
                    abs_mz_low=self.min_mz,
                    n_concatenated=self.n_concatenated,
                    concat_dt_idxs=concat_dt_idxs,
                    total_mass_window=self.total_mass_window,
                ))
            pmem(str(n_itr) + " End Factor " + str(i))
            n_itr += 1
        pmem(str(n_itr) + " Factor Initialization End")
        n_itr += 1
        self.factors = factors
        pmem(str(n_itr) + " Script End")
        # t = time.time()
        # print('Done: T+'+str(t-t0))


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
        factor_idx,
        n_factors,
        lows,
        highs,
        abs_mz_low,
        n_concatenated,
        concat_dt_idxs,
        total_mass_window,
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
        self.auc = sum(mz_data)
        self.factor_idx = factor_idx
        self.n_factors = n_factors
        self.lows = lows
        self.highs = highs
        self.abs_mz_low = abs_mz_low
        self.n_concatenated = n_concatenated
        self.concat_dt_idxs = concat_dt_idxs
        self.total_mass_window = total_mass_window

        ###Compute Instance Values###

        # integrate within expected peak bounds and create boolean mask of expected peak bounds called grate
        self.integrated_mz_data = []
        self.grate = np.resize(False, (len(self.mz_data)))
        for i, j in zip(self.lows, self.highs):
            self.integrated_mz_data.append(sum(self.mz_data[i:j]))
            self.grate[i:j] = True
        self.grate_sum = sum(self.mz_data[self.grate])

        self.max_rtdt = max(self.rts) * max(self.dts)
        self.outer_rtdt = sum(sum(np.outer(self.rts, self.dts)))

        # This can be a shared function
        self.integration_box_centers = []
        [
            self.integration_box_centers.append(i + ((j - i) / 2))
            for i, j in zip(self.lows, self.highs)
        ]

        self.integrated_mz_baseline = peakutils.baseline(
            np.asarray(self.integrated_mz_data),
            6)  # 6 degree curve seems to work well
        self.baseline_subtracted_integrated_mz = (self.integrated_mz_data -
                                                  self.integrated_mz_baseline)

        # this is a poor implementation, at least use list comprehensions TODO
        self.box_dist_avg = 0
        for i in range(1, len(self.integration_box_centers)):
            self.box_dist_avg += (self.integration_box_centers[i] -
                                  self.integration_box_centers[i - 1])
        self.box_dist_avg = self.box_dist_avg / (
            len(self.integration_box_centers) - 1)

        # Writes to self.isotope_clusters
        self.find_isotope_clusters(
            5, height=0.5
        )  # heuristic height value, should be high-level param TODO - Will require passage through DataTensor class

    # Uses find_window function to identify portions of the integrated mz dimension that look 'isotope-cluster-like', saves as Factor attribute
    def find_isotope_clusters(self, peak_width, **kwargs):

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

        self.isotope_clusters = []
        peaks, feature_dict = find_peaks(self.baseline_subtracted_integrated_mz,
                                         prominence=0.01,
                                         width=0.5)

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
                peaks, self.baseline_subtracted_integrated_mz)
            [ic_idxs.append(tup) for tup in height_filtered]
            cluster_idx = 0
            for tup in ic_idxs:
                integrated_indices = tup
                if integrated_indices != None:
                    #try:
                    newIC = IsotopeCluster(
                        charge_states=self.charge_states,
                        factor_mz_data=copy.deepcopy(self.mz_data),
                        source_file=self.source_file,
                        tensor_idx=self.tensor_idx,
                        timepoint_idx=self.timepoint_idx,
                        n_factors=self.n_factors,
                        factor_idx=self.factor_idx,
                        cluster_idx=cluster_idx,
                        low_idx=self.lows[integrated_indices[0]] -
                        math.ceil(self.box_dist_avg / 2),
                        high_idx=self.highs[integrated_indices[1]] +
                        math.ceil(self.box_dist_avg / 2),
                        lows=self.lows,
                        highs=self.highs,
                        grate=self.grate,
                        rts=self.rts,
                        dts=self.dts,
                        max_rtdt=self.max_rtdt,
                        outer_rtdt=self.outer_rtdt,
                        box_dist_avg=self.box_dist_avg,
                        abs_mz_low=self.abs_mz_low,
                        n_concatenated=self.n_concatenated,
                        concat_dt_idxs=self.concat_dt_idxs,
                        total_mass_window=self.total_mass_window,
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
        total_mass_window,
        lows,
        highs,
        grate,
        rts,
        dts,
        max_rtdt,
        outer_rtdt,
        box_dist_avg,
        abs_mz_low,
        n_concatenated,
        concat_dt_idxs,
    ):

        ###Set Attributes###

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
        self.lows = lows
        self.highs = highs
        self.grate = grate
        self.rts = rts
        self.dts = dts
        self.max_rtdt = max_rtdt
        self.outer_rtdt = outer_rtdt
        self.box_dist_avg = box_dist_avg
        self.abs_mz_low = abs_mz_low
        self.n_concatenated = n_concatenated
        self.concat_dt_idxs = concat_dt_idxs
        self.total_mass_window = total_mass_window

        ###Calculate Scoring Requirements###

        # Create array of expected peak positions
        self.integration_box_centers = []
        for i, j in zip(self.lows, self.highs):
            self.integration_box_centers.append(i + ((j - i) / 2))

        # prune factor_mz to get window around cluster that is consistent between charge-states
        self.cluster_mz_data = copy.deepcopy(self.factor_mz_data)
        self.cluster_mz_data[0:self.low_idx] = 0
        self.cluster_mz_data[self.high_idx:] = 0

        # integrate area of IC
        self.auc = sum(self.cluster_mz_data) * self.outer_rtdt

        # Values of isotope cluster that fall within the grate of expected peak bounds
        self.box_intensities = self.cluster_mz_data[self.grate]
        self.grate_sum = sum(self.box_intensities)

        # identify peaks and find error from expected peak positions using raw mz

        self.mz_peaks = find_peaks(self.factor_mz_data,
                                   distance=self.box_dist_avg)[0]
        self.max_peak_height = max(self.box_intensities)

        self.peak_error, self.peaks_chosen = self.find_peak_error(
            self.cluster_mz_data,
            self.mz_peaks,
            self.
            integration_box_centers[np.searchsorted(self.lows, self.low_idx):np.
                                    searchsorted(self.highs, self.high_idx)],
            self.max_peak_height,
        )

        # subtract baseline from IC mz values, recompute intrinsic values with new array
        self.baseline = peakutils.baseline(
            self.cluster_mz_data[self.low_idx:self.high_idx],
            6)  # 6 degree curve seems to work well
        self.baseline_subtracted_mz = self.cluster_mz_data
        self.baseline_subtracted_mz[self.low_idx:self.high_idx] = (
            self.cluster_mz_data[self.low_idx:self.high_idx] - self.baseline)
        self.baseline_auc = sum(self.baseline_subtracted_mz) * self.outer_rtdt
        self.log_baseline_auc = np.log(self.baseline_auc)
        self.baseline_box_intensities = self.baseline_subtracted_mz[self.grate]
        self.baseline_grate_sum = sum(self.baseline_box_intensities)
        self.baseline_max_peak_height = max(self.baseline_box_intensities)
        self.baseline_peak_error, self.baseline_peaks_chosen = self.find_peak_error(
            self.baseline_subtracted_mz,
            self.mz_peaks,
            self.
            integration_box_centers[np.searchsorted(self.lows, self.low_idx):np.
                                    searchsorted(self.highs, self.high_idx)],
            self.baseline_max_peak_height,
        )

        # create integrated mz array, indexed by integration box
        baseline_int_mz = []
        for lo, hi in zip(self.lows, self.highs):
            baseline_int_mz.append(sum(self.baseline_subtracted_mz[lo:hi]))
        self.baseline_integrated_mz = np.asarray(baseline_int_mz)

        # Cache int_mz and rt scoring values
        self.baseline_integrated_mz_norm = self.baseline_integrated_mz / np.linalg.norm(
            self.baseline_integrated_mz)
        self.baseline_integrated_mz_com = center_of_mass(
            self.baseline_integrated_mz)[
                0]  # COM in IC integrated bin dimension
        self.baseline_integrated_mz_std = (np.average(
            (np.arange(len(self.baseline_integrated_mz)) -
             self.baseline_integrated_mz_com)**2,
            weights=self.baseline_integrated_mz,
        )**0.5)

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
            self.abs_mz_com = self.find_mz_com(self.total_mass_window)
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
            self.baseline_grate_sum,  # Baseline-subtracted grate area sum (BGS)
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

    # uses internal
    def find_mz_com(self, tensor_mass_range):
        factor_mz_range = tensor_mass_range / self.charge_states[0]
        factor_mz_bin_step = factor_mz_range / len(self.factor_mz_data)
        left_mz_dist = (
            self.highs[math.floor(self.baseline_integrated_mz_com)] -
            self.integration_box_centers[math.floor(
                self.baseline_integrated_mz_com)]
        ) * factor_mz_bin_step  # MZ dist from center to right bound, times bin_to_mz factor
        right_mz_dist = (
            self.integration_box_centers[
                math.floor(self.baseline_integrated_mz_com) + 1] -
            self.lows[math.floor(self.baseline_integrated_mz_com) + 1]
        ) * factor_mz_bin_step  # MZ dist from left bound to center, times bin_to_mz factor

        if (self.baseline_integrated_mz_com -
                math.floor(self.baseline_integrated_mz_com)) <= 0.5:
            major = self.abs_mz_low + (factor_mz_bin_step *
                                       self.integration_box_centers[math.floor(
                                           self.baseline_integrated_mz_com)])
            minor = left_mz_dist * (self.baseline_integrated_mz_com -
                                    math.floor(self.baseline_integrated_mz_com))
            abs_mz_com = major + minor
        else:
            major = self.abs_mz_low + (
                factor_mz_bin_step *
                self.lows[math.floor(self.baseline_integrated_mz_com)] + 1)
            minor = right_mz_dist * (
                self.baseline_integrated_mz_com -
                math.floor(self.baseline_integrated_mz_com) - 0.5)
            abs_mz_com = major + minor

        return abs_mz_com

    # calculates sum of distances between nearest prominent peaks and expected peak centers in IC
    def find_peak_error(self, source, mz_peaks, integration_box_centers,
                        max_peak_height):

        peak_error = 0
        peaks_chosen = []
        peaks_total_height = 0
        match_idx = np.searchsorted(mz_peaks, integration_box_centers)
        if len(mz_peaks) > 0:
            for i in range(len(match_idx)):
                # handle peaks list of length 1
                if len(mz_peaks) == 1:
                    peak_error += abs(integration_box_centers[i] -
                                      mz_peaks[0]) * (source[mz_peaks[0]])
                    peaks_chosen.append(mz_peaks[0])
                    peaks_total_height += source[mz_peaks[0]]
                else:
                    # check if place to be inserted is leftmost of peaks
                    if match_idx[i] == 0:
                        peak_error += abs(integration_box_centers[i] -
                                          mz_peaks[match_idx[i]]) * (
                                              source[mz_peaks[match_idx[i]]])
                        peaks_chosen.append(mz_peaks[match_idx[i]])
                        peaks_total_height += source[mz_peaks[match_idx[i]]]
                    else:
                        # check if insertion position is rightmost of peaks
                        if match_idx[i] == len(mz_peaks):
                            peak_error += abs(integration_box_centers[i] -
                                              mz_peaks[-1]) * (
                                                  source[mz_peaks[-1]])
                            peaks_chosen.append(mz_peaks[-1])
                            peaks_total_height += source[mz_peaks[-1]]
                        else:
                            # handle case where distances between peaks are the same, pick biggest peak
                            if abs(integration_box_centers[i] -
                                   mz_peaks[match_idx[i]]) == abs(
                                       integration_box_centers[i] -
                                       mz_peaks[match_idx[i] - 1]):
                                peak_error += max([
                                    abs(integration_box_centers[i] -
                                        mz_peaks[match_idx[i]]) *
                                    (source[mz_peaks[match_idx[i]]]),
                                    abs(integration_box_centers[i] -
                                        mz_peaks[match_idx[i] - 1]) *
                                    (source[mz_peaks[match_idx[i] - 1]]),
                                ])
                                if abs(integration_box_centers[i] -
                                       mz_peaks[match_idx[i]]) * (
                                           source[mz_peaks[match_idx[i]]]
                                       ) > abs(integration_box_centers[i] -
                                               mz_peaks[match_idx[i] - 1]) * (
                                                   source[mz_peaks[match_idx[i]
                                                                   - 1]]):
                                    peaks_chosen.append(mz_peaks[match_idx[i]])
                                    peaks_total_height += source[mz_peaks[
                                        match_idx[i]]]
                                else:
                                    peaks_chosen.append(mz_peaks[match_idx[i] -
                                                                 1])
                                    peaks_total_height += source[mz_peaks[
                                        match_idx[i] - 1]]
                            else:
                                # only need to check left hand side differences because of left-hand default of searchsorted algorithm
                                # now check which peak is closer, left or right. This poses problems as there may be very close peaks which are not
                                # actually significant in height but which pass filtering.
                                if abs(integration_box_centers[i] -
                                       mz_peaks[match_idx[i]]) < abs(
                                           integration_box_centers[i] -
                                           mz_peaks[match_idx[i] - 1]):
                                    peak_error += abs(
                                        integration_box_centers[i] -
                                        mz_peaks[match_idx[i]]) * (
                                            source[mz_peaks[match_idx[i]]])
                                    peaks_chosen.append(mz_peaks[match_idx[i]])
                                    peaks_total_height += source[mz_peaks[
                                        match_idx[i]]]
                                else:
                                    peak_error += abs(
                                        integration_box_centers[i] -
                                        mz_peaks[match_idx[i] - 1]) * (
                                            source[mz_peaks[match_idx[i] - 1]])
                                    peaks_chosen.append(mz_peaks[match_idx[i] -
                                                                 1])
                                    peaks_total_height += source[mz_peaks[
                                        match_idx[i] - 1]]

            box_dist_total = 0
            for i in range(1, len(integration_box_centers)):
                box_dist_total += (integration_box_centers[i] -
                                   integration_box_centers[i - 1])

            if len(integration_box_centers) == 1:
                peak_error = (peak_error / peaks_total_height /
                              (box_dist_total / (len(integration_box_centers))))
            else:
                if len(integration_box_centers) > 1:
                    peak_error = (peak_error / peaks_total_height /
                                  (box_dist_total /
                                   (len(integration_box_centers) - 1)))
                else:
                    peak_error = 100000
            return peak_error, peaks_chosen
