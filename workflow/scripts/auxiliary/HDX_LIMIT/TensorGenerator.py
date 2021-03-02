import DataTensor

class TensorGenerator:
    # REWORK OF SerialTG FOR SINGLE TENSORS

    ###Class Attributes###
    hd_mass_diff = 1.006277
    c13_mass_diff = 1.00335

    def __init__(self, filename, timepoint_index, library_info, **kwargs):

        ###Set Instance Attributes###

        self.filename = filename
        self.timepoint_index = timepoint_index
        self.library_info = library_info

        if (
            kwargs is not None
        ):  # TODO: Default control values should be placed in a check statement below (if hasattr(self, 'name of value'): Do) do this for other kwargs fxns
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

        self.tensor = limit_read(self.filename)
        self.lib_idx = int(
            filename.split("/")[-1].split("_")[0]
        )  # expects format: path/to/{LibraryIdx}_{protName}_{tp}.cpickle.zlib
        self.name = self.library_info.iloc[self.lib_idx]["name"]
        self.max_peak_center = len(
            self.library_info.loc[self.library_info["name"] == self.name][
                "sequence"
            ].values[0]
        )
        self.total_isotopes = self.max_peak_center + self.high_mass_margin
        self.total_mass_window = self.low_mass_margin + self.total_isotopes

        self.est_peak_gaps = (
            [0]
            + list(np.linspace(self.c13_mass_diff, self.hd_mass_diff, 7))
            + [self.hd_mass_diff for x in range(self.total_isotopes - 8)]
        )
        self.cum_peak_gaps = np.cumsum(self.est_peak_gaps)

        i = self.lib_idx
        self.mz_centers = self.library_info["obs_mz"].values[i] + (
            self.cum_peak_gaps / self.library_info["charge"].values[i]
        )
        self.mz_lows = self.library_info["obs_mz"].values[i] - (
            self.low_mass_margin / self.library_info["charge"].values[i]
        )
        self.mz_highs = self.library_info["obs_mz"].values[i] + (
            self.total_isotopes / self.library_info["charge"].values[i]
        )

        self.low_lims = self.mz_centers * ((1000000.0 - self.ppm_radius) / 1000000.0)
        self.high_lims = self.mz_centers * ((1000000.0 + self.ppm_radius) / 1000000.0)

        # Instantitate DataTensor
        self.DataTensor = DataTensor(
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
        )

        self.DataTensor.lows = np.searchsorted(self.DataTensor.mz_labels, self.low_lims)
        self.DataTensor.highs = np.searchsorted(
            self.DataTensor.mz_labels, self.high_lims
        )
        # Consider separating factorize from init
        # self.DataTensor.factorize(gauss_params=(3,1))