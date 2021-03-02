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
        self.cluster_mz_data[0 : self.low_idx] = 0
        self.cluster_mz_data[self.high_idx :] = 0

        # integrate area of IC
        self.auc = sum(self.cluster_mz_data) * self.outer_rtdt

        # Values of isotope cluster that fall within the grate of expected peak bounds
        self.box_intensities = self.cluster_mz_data[self.grate]
        self.grate_sum = sum(self.box_intensities)

        # identify peaks and find error from expected peak positions using raw mz

        self.mz_peaks = sp.signal.find_peaks(
            self.factor_mz_data, distance=self.box_dist_avg
        )[0]
        self.max_peak_height = max(self.box_intensities)

        self.peak_error, self.peaks_chosen = self.find_peak_error(
            self.cluster_mz_data,
            self.mz_peaks,
            self.integration_box_centers[
                np.searchsorted(self.lows, self.low_idx) : np.searchsorted(
                    self.highs, self.high_idx
                )
            ],
            self.max_peak_height,
        )

        # subtract baseline from IC mz values, recompute intrinsic values with new array
        self.baseline = peakutils.baseline(
            self.cluster_mz_data[self.low_idx : self.high_idx], 6
        )  # 6 degree curve seems to work well
        self.baseline_subtracted_mz = self.cluster_mz_data
        self.baseline_subtracted_mz[self.low_idx : self.high_idx] = (
            self.cluster_mz_data[self.low_idx : self.high_idx] - self.baseline
        )
        self.baseline_auc = sum(self.baseline_subtracted_mz) * self.outer_rtdt
        self.log_baseline_auc = np.log(self.baseline_auc)
        self.baseline_box_intensities = self.baseline_subtracted_mz[self.grate]
        self.baseline_grate_sum = sum(self.baseline_box_intensities)
        self.baseline_max_peak_height = max(self.baseline_box_intensities)
        self.baseline_peak_error, self.baseline_peaks_chosen = self.find_peak_error(
            self.baseline_subtracted_mz,
            self.mz_peaks,
            self.integration_box_centers[
                np.searchsorted(self.lows, self.low_idx) : np.searchsorted(
                    self.highs, self.high_idx
                )
            ],
            self.baseline_max_peak_height,
        )

        # create integrated mz array, indexed by integration box
        baseline_int_mz = []
        for lo, hi in zip(self.lows, self.highs):
            baseline_int_mz.append(sum(self.baseline_subtracted_mz[lo:hi]))
        self.baseline_integrated_mz = np.asarray(baseline_int_mz)

        # Cache int_mz and rt scoring values
        self.baseline_integrated_mz_norm = self.baseline_integrated_mz / np.linalg.norm(
            self.baseline_integrated_mz
        )
        self.baseline_integrated_mz_com = sp.ndimage.measurements.center_of_mass(
            self.baseline_integrated_mz
        )[
            0
        ]  # COM in IC integrated bin dimension
        self.baseline_integrated_mz_std = (
            np.average(
                (
                    np.arange(len(self.baseline_integrated_mz))
                    - self.baseline_integrated_mz_com
                )
                ** 2,
                weights=self.baseline_integrated_mz,
            )
            ** 0.5
        )

        self.rt_norm = self.rts / np.linalg.norm(self.rts)
        self.rt_com = sp.ndimage.measurements.center_of_mass(self.rts)[0]

        # Cache DT values
        # If DT is concatenated, return list of coms and norms of single rts relative to bin numbers, a single_dt distribution starts at 0. If only one charge state, return list of len=1
        if self.concat_dt_idxs is not None:
            single_dts = []
            # generate list of single dts
            single_dts.append(self.dts[: self.concat_dt_idxs[0]])
            for i in range(len(self.charge_states) - 1):
                single_dts.append(
                    self.dts[self.concat_dt_idxs[i] : self.concat_dt_idxs[i + 1]]
                )

            self.single_dts = single_dts
            self.dt_coms = [
                sp.ndimage.measurements.center_of_mass(dt)[0] for dt in single_dts
            ]
            self.dt_norms = [dt / np.linalg.norm(dt) for dt in single_dts]
        else:
            self.dt_coms = [sp.ndimage.measurements.center_of_mass(self.dts)[0]]
            self.dt_norms = [self.dts / np.linalg.norm(self.dts)]

        if self.n_concatenated == 1:
            self.abs_mz_com = self.find_mz_com(self.total_mass_window)
        else:
            self.abs_mz_com = "Concatenated, N/A, see IC.baseline_integrated_mz_com"

        # format useful values to be read by pandas
        self.info_tuple = (
            self.source_file,  # Filename of data used to create parent DataTensor
            self.tensor_idx,  # Library master list row of parent-DataTensor
            self.n_factors,  # Number of factors in parent decomposition
            self.factor_idx,  # Index of IC parent-factor in DataTensor.factors[]
            self.cluster_idx,  # Index of IC in parent-factor.isotope_clusters[]
            self.charge_states,  # List of charge states in IC
            self.n_concatenated,  # number of source tensors IC parent-DataTensor was made from
            self.low_idx,  # Low bin index corresponding to Factor-level bins
            self.high_idx,  # High bin index corresponding to Factor-level bins
            self.baseline_auc,  # Baseline-subtracted AUC (BAUC)
            self.baseline_grate_sum,  # Baseline-subtracted grate area sum (BGS)
            self.baseline_peak_error,  # Baseline-subtracted version of peak-error (BPE)
            self.baseline_integrated_mz_com,  # Center of mass in added mass units
            self.abs_mz_com,  # IsotopicCluster center-of-mass in absolute m/z dimension
            self.rts,  # Array of RT values
            self.dts,  # Array of DT values, if a tensor is concatenated,this is taken from the last tensor in the list, can be seen in tensor_idx
            np.arange(0, len(self.baseline_integrated_mz), 1),
            self.baseline_integrated_mz,  # Array of baseline-subtracted integrated mz intensity values
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
            self.highs[math.floor(self.baseline_integrated_mz_com)]
            - self.integration_box_centers[math.floor(self.baseline_integrated_mz_com)]
        ) * factor_mz_bin_step  # MZ dist from center to right bound, times bin_to_mz factor
        right_mz_dist = (
            self.integration_box_centers[
                math.floor(self.baseline_integrated_mz_com) + 1
            ]
            - self.lows[math.floor(self.baseline_integrated_mz_com) + 1]
        ) * factor_mz_bin_step  # MZ dist from left bound to center, times bin_to_mz factor

        if (
            self.baseline_integrated_mz_com
            - math.floor(self.baseline_integrated_mz_com)
        ) <= 0.5:
            major = self.abs_mz_low + (
                factor_mz_bin_step
                * self.integration_box_centers[
                    math.floor(self.baseline_integrated_mz_com)
                ]
            )
            minor = left_mz_dist * (
                self.baseline_integrated_mz_com
                - math.floor(self.baseline_integrated_mz_com)
            )
            abs_mz_com = major + minor
        else:
            major = self.abs_mz_low + (
                factor_mz_bin_step
                * self.lows[math.floor(self.baseline_integrated_mz_com)]
                + 1
            )
            minor = right_mz_dist * (
                self.baseline_integrated_mz_com
                - math.floor(self.baseline_integrated_mz_com)
                - 0.5
            )
            abs_mz_com = major + minor

        return abs_mz_com

    # calculates sum of distances between nearest prominent peaks and expected peak centers in IC
    def find_peak_error(
        self, source, mz_peaks, integration_box_centers, max_peak_height
    ):

        peak_error = 0
        peaks_chosen = []
        peaks_total_height = 0
        match_idx = np.searchsorted(mz_peaks, integration_box_centers)
        if len(mz_peaks) > 0:
            for i in range(len(match_idx)):
                # handle peaks list of length 1
                if len(mz_peaks) == 1:
                    peak_error += abs(integration_box_centers[i] - mz_peaks[0]) * (
                        source[mz_peaks[0]]
                    )
                    peaks_chosen.append(mz_peaks[0])
                    peaks_total_height += source[mz_peaks[0]]
                else:
                    # check if place to be inserted is leftmost of peaks
                    if match_idx[i] == 0:
                        peak_error += abs(
                            integration_box_centers[i] - mz_peaks[match_idx[i]]
                        ) * (source[mz_peaks[match_idx[i]]])
                        peaks_chosen.append(mz_peaks[match_idx[i]])
                        peaks_total_height += source[mz_peaks[match_idx[i]]]
                    else:
                        # check if insertion position is rightmost of peaks
                        if match_idx[i] == len(mz_peaks):
                            peak_error += abs(
                                integration_box_centers[i] - mz_peaks[-1]
                            ) * (source[mz_peaks[-1]])
                            peaks_chosen.append(mz_peaks[-1])
                            peaks_total_height += source[mz_peaks[-1]]
                        else:
                            # handle case where distances between peaks are the same, pick biggest peak
                            if abs(
                                integration_box_centers[i] - mz_peaks[match_idx[i]]
                            ) == abs(
                                integration_box_centers[i] - mz_peaks[match_idx[i] - 1]
                            ):
                                peak_error += max(
                                    [
                                        abs(
                                            integration_box_centers[i]
                                            - mz_peaks[match_idx[i]]
                                        )
                                        * (source[mz_peaks[match_idx[i]]]),
                                        abs(
                                            integration_box_centers[i]
                                            - mz_peaks[match_idx[i] - 1]
                                        )
                                        * (source[mz_peaks[match_idx[i] - 1]]),
                                    ]
                                )
                                if abs(
                                    integration_box_centers[i] - mz_peaks[match_idx[i]]
                                ) * (source[mz_peaks[match_idx[i]]]) > abs(
                                    integration_box_centers[i]
                                    - mz_peaks[match_idx[i] - 1]
                                ) * (
                                    source[mz_peaks[match_idx[i] - 1]]
                                ):
                                    peaks_chosen.append(mz_peaks[match_idx[i]])
                                    peaks_total_height += source[mz_peaks[match_idx[i]]]
                                else:
                                    peaks_chosen.append(mz_peaks[match_idx[i] - 1])
                                    peaks_total_height += source[
                                        mz_peaks[match_idx[i] - 1]
                                    ]
                            else:
                                # only need to check left hand side differences because of left-hand default of searchsorted algorithm
                                # now check which peak is closer, left or right. This poses problems as there may be very close peaks which are not
                                # actually significant in height but which pass filtering.
                                if abs(
                                    integration_box_centers[i] - mz_peaks[match_idx[i]]
                                ) < abs(
                                    integration_box_centers[i]
                                    - mz_peaks[match_idx[i] - 1]
                                ):
                                    peak_error += abs(
                                        integration_box_centers[i]
                                        - mz_peaks[match_idx[i]]
                                    ) * (source[mz_peaks[match_idx[i]]])
                                    peaks_chosen.append(mz_peaks[match_idx[i]])
                                    peaks_total_height += source[mz_peaks[match_idx[i]]]
                                else:
                                    peak_error += abs(
                                        integration_box_centers[i]
                                        - mz_peaks[match_idx[i] - 1]
                                    ) * (source[mz_peaks[match_idx[i] - 1]])
                                    peaks_chosen.append(mz_peaks[match_idx[i] - 1])
                                    peaks_total_height += source[
                                        mz_peaks[match_idx[i] - 1]
                                    ]

            box_dist_total = 0
            for i in range(1, len(integration_box_centers)):
                box_dist_total += (
                    integration_box_centers[i] - integration_box_centers[i - 1]
                )

            if len(integration_box_centers) == 1:
                peak_error = (
                    peak_error
                    / peaks_total_height
                    / (box_dist_total / (len(integration_box_centers)))
                )
            else:
                if len(integration_box_centers) > 1:
                    peak_error = (
                        peak_error
                        / peaks_total_height
                        / (box_dist_total / (len(integration_box_centers) - 1))
                    )
                else:
                    peak_error = 100000
            return peak_error, peaks_chosen
