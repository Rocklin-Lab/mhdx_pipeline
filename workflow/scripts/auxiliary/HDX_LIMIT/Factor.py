import IsotopeCluster

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
            np.asarray(self.integrated_mz_data), 6
        )  # 6 degree curve seems to work well
        self.baseline_subtracted_integrated_mz = (
            self.integrated_mz_data - self.integrated_mz_baseline
        )

        # this is a poor implementation, at least use list comprehensions TODO
        self.box_dist_avg = 0
        for i in range(1, len(self.integration_box_centers)):
            self.box_dist_avg += (
                self.integration_box_centers[i] - self.integration_box_centers[i - 1]
            )
        self.box_dist_avg = self.box_dist_avg / (len(self.integration_box_centers) - 1)

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
        peaks, feature_dict = sp.signal.find_peaks(
            self.baseline_subtracted_integrated_mz, prominence=0.01, width=0.5
        )

        if len(peaks) == 0:
            return
        else:
            ic_idxs = [
                (feature_dict["left_bases"][i], feature_dict["right_bases"][i])
                for i in range(len(peaks))
                if feature_dict["left_bases"][i] < feature_dict["right_bases"][i]
                if feature_dict["right_bases"][i] - feature_dict["left_bases"][i] > 4
            ]
            # ic_idxs = [(feature_dict['left_bases'][i], feature_dict['left_bases'][i+1]) if feature_dict['left_bases'][i] < feature_dict['left_bases'][i+1] else (feature_dict['left_bases'][i], feature_dict['left_bases'][i]+6) for i in range(len(out[0])-1)]
            if len(peaks) > 1:
                ic_idxs.append(
                    (feature_dict["left_bases"][0], feature_dict["right_bases"][-1])
                )  # Create default ic from first left base to last right base
            height_filtered = rel_height_peak_bounds(
                peaks, self.baseline_subtracted_integrated_mz
            )
            [ic_idxs.append(tup) for tup in height_filtered]
            cluster_idx = 0
            for tup in ic_idxs:
                integrated_indices = tup
                if integrated_indices != None:
                    try:
                        newIC = IsotopeCluster(
                            charge_states=self.charge_states,
                            factor_mz_data=copy.deepcopy(self.mz_data),
                            source_file=self.source_file,
                            tensor_idx=self.tensor_idx,
                            timepoint_idx=self.timepoint_idx,
                            n_factors=self.n_factors,
                            factor_idx=self.factor_idx,
                            cluster_idx=cluster_idx,
                            low_idx=self.lows[integrated_indices[0]]
                            - math.ceil(self.box_dist_avg / 2),
                            high_idx=self.highs[integrated_indices[1]]
                            + math.ceil(self.box_dist_avg / 2),
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
                        if (
                            newIC.baseline_peak_error / newIC.baseline_auc < 0.2
                        ):  # TODO: HARDCODE
                            self.isotope_clusters.append(newIC)
                            cluster_idx += 1
                    except:
                        print("ic index out of bounds: " + str(integrated_indices))
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
            if (
                array[idx] < array[peak_idx] / 5
            ):  # Peak is likely not an IC if peak > 5 x neighbors
                win_high = idx
                rflag = False

        while rflag:
            # make sure looking ahead won't throw error
            if idx + 1 < len(array):
                # if idx+1 goes down, and its height is greater than 20% of the max peak
                if array[idx + 1] < array[idx] and array[idx + 1] > array[peak_idx] / 5:
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
                    array[idx - 1] < array[idx] and array[idx - 1] > array[peak_idx] / 5
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

