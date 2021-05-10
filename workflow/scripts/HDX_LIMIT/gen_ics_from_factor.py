import numpy as np
from scipy.signal import find_peaks

# def find_isotope_clusters(self, peak_width, **kwargs):
#     def rel_height_peak_bounds(centers, int_mz, bound=20):
#         out = []
#         baseline = max(int_mz) * 0.15  # TODO: HARDCODE
#         for center in centers:
#             if int_mz[center] > baseline:
#                 i, j = center, center
#                 cutoff = int_mz[center] * (bound / 100)
#                 while center - i <= 10 and i - 1 != -1:
#                     i -= 1
#                     if int_mz[i] < cutoff:
#                         break
#                 while j - center <= 10 and j + 1 != len(int_mz):
#                     j += 1
#                     if int_mz[j] < cutoff:
#                         break
#                 out.append((i, j))
#         return out
#
#     self.isotope_clusters = []
#     peaks, feature_dict = find_peaks(self.baseline_subtracted_integrated_mz,
#                                      prominence=0.01,
#                                      width=0.5)
#
#     if len(peaks) == 0:
#         return
#     else:
#         ic_idxs = [
#             (feature_dict["left_bases"][i], feature_dict["right_bases"][i])
#             for i in range(len(peaks))
#             if
#             feature_dict["left_bases"][i] < feature_dict["right_bases"][i]
#             if feature_dict["right_bases"][i] -
#             feature_dict["left_bases"][i] > 4
#         ]
#         # ic_idxs = [(feature_dict['left_bases'][i], feature_dict['left_bases'][i+1]) if feature_dict['left_bases'][i] < feature_dict['left_bases'][i+1] else (feature_dict['left_bases'][i], feature_dict['left_bases'][i]+6) for i in range(len(out[0])-1)]
#         if len(peaks) > 1:
#             ic_idxs.append(
#                 (feature_dict["left_bases"][0],
#                  feature_dict["right_bases"][-1])
#             )  # Create default ic from first left base to last right base
#         height_filtered = rel_height_peak_bounds(
#             peaks, self.baseline_subtracted_integrated_mz)
#         [ic_idxs.append(tup) for tup in height_filtered]
#         cluster_idx = 0
#         for integrated_indices in ic_idxs:
#             if integrated_indices != None:
#                 # try:
#
#                 newIC = IsotopeCluster(
#                     charge_states=self.charge_states,
#                     factor_mz_data=copy.deepcopy(self.mz_data),
#                     source_file=self.source_file,
#                     tensor_idx=self.tensor_idx,
#                     timepoint_idx=self.timepoint_idx,
#                     n_factors=self.n_factors,
#                     factor_idx=self.factor_idx,
#                     cluster_idx=cluster_idx,
#                     low_idx=self.bins_per_isotope_peak * integrated_indices[0],
#                     high_idx=self.bins_per_isotope_peak * (integrated_indices[1] + 1),
#                     rts=self.rts,
#                     dts=self.dts,
#                     retention_labels=self.retention_labels,
#                     drift_labels=self.drift_labels,
#                     mz_labels=self.mz_labels,
#                     bins_per_isotope_peak=self.bins_per_isotope_peak,
#                     max_rtdt=self.max_rtdt,
#                     outer_rtdt=self.outer_rtdt,
#                     n_concatenated=self.n_concatenated,
#                     concat_dt_idxs=self.concat_dt_idxs,
#                 )
#                 if (newIC.baseline_peak_error / newIC.baseline_auc <
#                     0.2):  # TODO: HARDCODE
#                     self.isotope_clusters.append(newIC)
#                     cluster_idx += 1
#                 # except:
#                 # print("ic index out of bounds: " + str(integrated_indices))
#         return


def find_isotope_cluster(factors):
    """
    find isotope cluster from factor data
    :param factors: factors dictionary from factor data dictionary
    :return:
    """

    # use integrated mz array to find peaks

    for factor in factors:
        integrated_mz = factor['factor_integrated_mz']
