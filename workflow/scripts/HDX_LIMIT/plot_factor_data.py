import seaborn as sns



def plot_factor_row(fig, gs, retention_labels, drift_lables, mz_labels, bins_per_isotope_peak, row_number, factor=None,
                    tensor3=None, name=''):
    """
    plot the factor data (3d rt, dt, mz), mz data, integrated mz data
    :param retention_labels: retention labels
    :param drift_lables: drift labels
    :param mz_labels: mz labels
    :param bins_per_isotope_peak: number of bins per isotope peak
    :param row_number: row number for plotting
    :param factor: factor
    :param tensor3: 3d grid (rt, dt, mz)
    :param name: name of the plot
    :return:
    """

    if factor != None:
        factor_rt_dt_grid = np.multiply.outer(factor['factor_dt'], factor['factor_rt'])
        factor_mz = factor['factor_mz']
        factor_integrated_mz = factor['factor_integrated_mz']
    else:
        factor_rt_dt_grid = np.sum(tensor3, axis=2).T
        factor_mz = np.sum(tensor3, axis=(0, 1))
        factor_integrated_mz = np.sum(np.reshape(factor_mz, (-1, bins_per_isotope_peak)), axis=1)

    ax = fig.add_subplot(gs[row_number, 0])

    sns.heatmap(factor_rt_dt_grid, cbar=False, cmap='Blues')
    plt.xlabels()


def plot_factor_data(retention_labels, drift_labels, mz_labels, bins_per_isotope_peak, tensor3, factors,
                     gauss_filter_params=(3,1), title='', output_path=None):
    """
    plot factor data
    :param retention_labels: retention time labels
    :param drift_labels: drift time labels
    :param mz_labels: mz labels
    :param bins_per_isotope_peak: number of bins per isotope peak
    :param tensor3: 3d tensor (rt, dt, mz) grid
    :param factors: factor
    :param gauss_filter_params: gaussian filter params in tuple(rt sigma, dt sigma)
    :param title: title for the plot
    :param output_path: output path
    :return: None
    """
