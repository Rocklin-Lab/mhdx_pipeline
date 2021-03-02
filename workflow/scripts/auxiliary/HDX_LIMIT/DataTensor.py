import Factor

class DataTensor:
    def __init__(
        self,
        source_file,
        tensor_idx,
        timepoint_idx,
        name,
        total_mass_window,
        n_concatenated,
        charge_states,
        **kwargs
    ):

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
                    self.int_seq_out_float, (len(self.rts), len(self.dts), 50)
                )
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
            if not all(
                "dts" in kwargs
                and "rts" in kwargs
                and "lows" in kwargs
                and "highs" in kwargs
                and "concatenated_grid" in kwargs
                and "abs_mz_low" in kwargs
                and "concat_dt_idxs" in kwargs
            ):

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
            (len(retention_labels), len(drift_labels), len(mz_labels))
        )

        scan = 0
        for i in range(len(retention_labels)):
            for j in range(len(drift_labels)):
                mz_indices = np.searchsorted(mz_labels, sparse_data[scan][:, 0])
                large_array = np.zeros((len(sparse_data[scan]), len(mz_labels)))
                large_array[
                    np.arange(len(sparse_data[scan])), mz_indices
                ] = sparse_data[scan][:, 1]
                tensor3_out[i][j] = np.sum(large_array, axis=0)
                scan += 1

        return (retention_labels, drift_labels, mz_labels, min_mz, max_mz, tensor3_out)

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
                for _ in range(int(round(floor * ip)), int(round(ceil * ip)) + 1):
                    if _ % int(5) == 0:
                        a = float(_) / float(ip)
                        y = i * math.exp(-1 * ((mz - a) * (mz - a)) / (2 * s2))
                        tmp[a] += y
            out.append(np.asarray([[key, tmp[key]] for key in list(tmp.keys())]))
        return np.asarray(out)

    # Takes tensor input and gaussian filter parameters, outputs filtered data
    def gauss(self, grid, rt_sig=3, dt_sig=1):

        gauss_grid = np.zeros(np.shape(grid))
        for i in range(np.shape(grid)[2]):
            gauss_grid[:, :, i] = gaussian_filter(grid[:, :, i], (rt_sig, dt_sig))
        return gauss_grid

    def interpolate(self, grid_in, new_mz_len, gauss_params=None):
        # Takes length of mz_bins to interpolated to, and optional gaussian filter parameters
        # Returns the interpolated tensor, length of interpolated axis, interpolated low_lims and high_lims

        if gauss_params != None:
            grid = self.gauss(grid_in, gauss_params[0], gauss_params[1])
        else:
            grid = grid_in

        test_points = []
        z_axis = np.clip(
            np.linspace(0, np.shape(grid)[2], new_mz_len), 0, np.shape(grid)[2] - 1
        )
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
            points=[x, y, z], values=grid
        )
        interpolated_out = interpolation_function(test_points)
        interpolated_out = np.reshape(
            interpolated_out, (np.shape(grid)[0], np.shape(grid)[1], new_mz_len)
        )

        interpolated_bin_mzs = np.linspace(
            self.mz_bin_low, self.mz_bin_high, new_mz_len
        )
        interpolated_low_lims = np.searchsorted(
            interpolated_bin_mzs, self.mz_labels[self.lows]
        )
        interpolated_high_lims = np.searchsorted(
            interpolated_bin_mzs, self.mz_labels[self.highs]
        )

        return [interpolated_out, interpolated_low_lims, interpolated_high_lims]

    def factorize(self, n_factors=13, new_mz_len=None, gauss_params=None):
        # Test factorization starting at n_factors = 15 and counting down, keep factorization that has no factors with correlation greater than 0.2 in any dimension.

        def corr_check(factors, cutoff):
            # Checks scipy non_negatve_parafac output factors for inter-factor (off-diagonal) correlations > cutoff, returns True if all values are < cutoff

            a = np.minimum(
                np.minimum(np.corrcoef(factors[0].T), np.corrcoef(factors[1].T)),
                np.corrcoef(factors[2].T),
            )

            if any(a[np.where(~np.eye(a.shape[0], dtype=bool))] > cutoff):
                return False
            else:
                return True

        def pmem(id_str):
            process = psutil.Process(os.getpid())
            print(
                id_str
                + " Process Memory (GB): "
                + str(process.memory_info().rss / 1024 / 1024 / 1024)
            )

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
                    grid = self.gauss(
                        self.full_grid_out, gauss_params[0], gauss_params[1]
                    )
                else:
                    grid = self.full_grid_out
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
                print(
                    "All n-factors failed for Index: "
                    + str(self.name)
                    + ", keeping 1 factor decomposition."
                )
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
                )
            )
            pmem(str(n_itr) + " End Factor " + str(i))
            n_itr += 1
        pmem(str(n_itr) + " Factor Initialization End")
        n_itr += 1
        self.factors = factors
        pmem(str(n_itr) + " Script End")
        # t = time.time()
        # print('Done: T+'+str(t-t0))