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
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter



def plot_factor_row(fig, gs, retention_labels, drift_labels, mz_labels, bins_per_isotope_peak, row_number, factor=None,
                    tensor3=None, name='', from_dict=True):
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
    :return: None
    """

    if factor != None:
        if from_dict:
            factor_rt_dt_grid = np.multiply.outer(factor['factor_dt'], factor['factor_rt'])
            factor_mz = factor['factor_mz']
            factor_integrated_mz = factor['factor_integrated_mz']
        else:
            factor_rt_dt_grid = np.multiply.outer(factor.dts, factor.rts)
            factor_mz = factor.mz_data
            factor_integrated_mz = factor.integrated_mz_data
    else:
        factor_rt_dt_grid = np.sum(tensor3, axis=2).T
        factor_mz = np.sum(tensor3, axis=(0, 1))
        factor_integrated_mz = np.sum(np.reshape(factor_mz, (-1, bins_per_isotope_peak)), axis=1)

    # plot rt dt heat map
    ax = fig.add_subplot(gs[row_number, 0])
    sns.heatmap(factor_rt_dt_grid, cbar=False, cmap='Blues')
    plt.xlabel('retention (minutes)', labelpad=-10)
    plt.ylabel('%s\n\ndrift (ms)' % name)
    ax.set_xticks([0, len(retention_labels)])
    ax.set_yticks([0, len(drift_labels)])
    ax.set_xticklabels(['        %.2f' % retention_labels[0], '%.2f        ' % retention_labels[-1]],
                       rotation='horizontal')
    ax.set_yticklabels(['%.1f' % drift_labels[0], '%.1f' % drift_labels[-1]], rotation='horizontal')
    ax.tick_params(length=3, pad=3)

    # plot mz data

    padded_mz_labels = []
    padded_factor_mz = []
    mz_label_spacing = mz_labels[1] - mz_labels[0]

    for i in range(0, len(mz_labels), bins_per_isotope_peak):
        padded_mz_labels.append(mz_labels[i] - mz_label_spacing)
        padded_mz_labels += list(mz_labels[i:i + bins_per_isotope_peak])
        padded_mz_labels.append(padded_mz_labels[-1] + mz_label_spacing)

        padded_factor_mz.append(0)
        padded_factor_mz += list(factor_mz[i:i + bins_per_isotope_peak])
        padded_factor_mz.append(0)

    ax = fig.add_subplot(gs[row_number, 1])
    plt.plot(padded_mz_labels, padded_factor_mz, linewidth=0.5)
    ax.set_yticks([])
    ax.tick_params(length=3, pad=3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # plot integrated mz data

    ax = fig.add_subplot(gs[row_number, 2])
    plt.plot(factor_integrated_mz, linewidth=0.5, marker='o', markersize=3)
    plt.xticks(range(0, len(factor_integrated_mz) + 5, 5))
    ax.set_yticks([])
    ax.tick_params(length=3, pad=3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.grid(axis='x', linewidth=0.25)


def plot_factor_data(retention_labels, drift_labels, mz_labels, bins_per_isotope_peak, tensor3, factors, output_path,
                     gauss_filter_params=(3, 1), title='', from_dict=True):
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

    if from_dict:
        print('plotting factor data from dictionary')
    else:
        print('plotting factor data from data tensor')

    n_factors = len(factors)
    if n_factors == 1:
        total_pdf_rows = n_factors + 2
    else:
        total_pdf_rows = n_factors + 8

    fig = plt.figure(figsize=(15, total_pdf_rows * 1.6))
    gs = gridspec.GridSpec(ncols=3, nrows=total_pdf_rows, figure=fig)

    # plot the raw data

    plot_factor_row(fig=fig,
                    gs=gs,
                    retention_labels=retention_labels,
                    drift_labels=drift_labels,
                    mz_labels=mz_labels,
                    bins_per_isotope_peak=bins_per_isotope_peak,
                    row_number=0,
                    factor=None,
                    tensor3=tensor3,
                    name='Raw')


    # plot the gaussian filtered raw data

    gauss_filtered_tensor3 = gaussian_filter(tensor3, sigma=[gauss_filter_params[0], gauss_filter_params[1], 0])

    plot_factor_row(fig=fig,
                    gs=gs,
                    retention_labels=retention_labels,
                    drift_labels=drift_labels,
                    mz_labels=mz_labels,
                    bins_per_isotope_peak=bins_per_isotope_peak,
                    row_number=1,
                    factor=None,
                    tensor3=gauss_filtered_tensor3,
                    name='Gaussian Filtered')


    # plot factor data

    for num, factor in enumerate(factors):

        if from_dict:

            plot_factor_row(fig=fig,
                            gs=gs,
                            retention_labels=retention_labels,
                            drift_labels=drift_labels,
                            mz_labels=mz_labels,
                            bins_per_isotope_peak=bins_per_isotope_peak,
                            row_number=num+2,
                            factor=factor,
                            tensor3=None,
                            name='Factor %s ' % factor['factor_num'],
                            from_dict=from_dict)
        else:
            plot_factor_row(fig=fig,
                            gs=gs,
                            retention_labels=retention_labels,
                            drift_labels=drift_labels,
                            mz_labels=mz_labels,
                            bins_per_isotope_peak=bins_per_isotope_peak,
                            row_number=num + 2,
                            factor=factor,
                            tensor3=None,
                            name='Factor %s ' % factor.factor_idx,
                            from_dict=from_dict)



    # additional plots if there are more than one factors present

    if n_factors > 1:

        # generate factor masses, factor rts, factor dts, and factor mzs
        total_factor_masses = []
        factor_dts = []
        factor_rts = []
        factor_mzs = []

        for num, factor in enumerate(factors):
            if from_dict:
                factor_dts.append(factor['factor_dt'] / max(factor['factor_dt']))
                factor_rts.append(factor['factor_rt'] / max(factor['factor_rt']))
                factor_mzs.append(factor['factor_mz'] / max(factor['factor_mz']))
                total_factor_masses.append(sum(factor['factor_integrated_mz']))
            else:
                factor_dts.append(factor.dts / max(factor.dts))
                factor_rts.append(factor.rts / max(factor.rts))
                factor_mzs.append(factor.mz_data / max(factor.mz_data))
                total_factor_masses.append(sum(factor.integrated_mz_data))


        total_factor_masses = np.array(total_factor_masses)
        factor_dts = np.array(factor_dts).T
        factor_rts = np.array(factor_rts).T
        factor_mzs = np.array(factor_mzs).T

        # plot factor rts
        ax = fig.add_subplot(gs[-6:-3, 0])
        sns.heatmap(factor_rts, cbar=False, cmap='Blues')
        plt.xlabel('Factor index', labelpad=1)
        plt.ylabel('Retention time (minutes)')
        ax.tick_params(length=3, pad=3)
        ytick_indices = np.searchsorted(np.arange(len(retention_labels)), np.linspace(0, len(retention_labels) - 1, 7))
        ax.set_yticks(ytick_indices)
        ax.set_yticklabels(['%.2f' % x for x in np.array(retention_labels)[ytick_indices]])

        # plot factor dts
        ax = fig.add_subplot(gs[-6:-3, 1])
        sns.heatmap(factor_dts, cbar=False, cmap='Blues')
        plt.xlabel('Factor index', labelpad=1)
        plt.ylabel('Drift time bin')
        plt.yticks([])
        ax.tick_params(length=3, pad=3)

        # plot overall mass of each factor
        ax = fig.add_subplot(gs[-6:-3, 2])
        plt.bar(range(n_factors), total_factor_masses / sum(total_factor_masses))
        plt.xlabel('Factor index', labelpad=1)
        plt.ylabel('Relative factor mass')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticklabels([])
        ax.set_xticks(range(n_factors))
        ax.set_xticklabels(range(n_factors))
        ax.tick_params(length=3, pad=3)

        annot_size = 12
        if n_factors > 8: annot_size = 10
        if n_factors > 12: annot_size = 8

        # plot rt correlations between factors
        ax = fig.add_subplot(gs[-3:, 0])
        sns.heatmap(np.corrcoef(factor_rts.T), vmin=0, vmax=1, annot=True, fmt='.2f', cmap='Blues', cbar=False,
                    annot_kws={"size": annot_size})
        ax.tick_params(length=3, pad=3)
        plt.xlabel('Retention correlation', labelpad=1)
        plt.ylabel('Factor index', labelpad=1)

        # plot dt correlations between factors
        ax = fig.add_subplot(gs[-3:, 1])
        sns.heatmap(np.corrcoef(factor_dts.T), vmin=0, vmax=1, annot=True, fmt='.2f', cmap='Blues', cbar=False,
                    annot_kws={"size": annot_size})
        ax.tick_params(length=3, pad=3)
        plt.xlabel('Drift time correlation', labelpad=1)

        # plot min correlations between rts, dts, and mzs
        ax = fig.add_subplot(gs[-3:, 2])
        sns.heatmap(
            np.minimum(np.minimum(np.corrcoef(factor_rts.T), np.corrcoef(factor_dts.T)), np.corrcoef(factor_mzs.T)),
            vmin=0, vmax=1, annot=True, fmt='.2f', cbar=False, cmap='Blues', annot_kws={"size": annot_size})
        ax.tick_params(length=3, pad=3)
        maxcorr = np.max(
            np.minimum(np.minimum(np.corrcoef(factor_rts.T), np.corrcoef(factor_dts.T)), np.corrcoef(factor_mzs.T))[
                np.triu_indices(n_factors, 1)])
        plt.xlabel('Min (RT, DT, mz intensity) correlation (max: %.2f)' % maxcorr, labelpad=1)


    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.1, top=0.95)
    plt.suptitle(title)

    plt.savefig(output_path)
    plt.close()


def plot_factor_data_from_data_tensor(data_tensor, output_path):

    title = ''
    if hasattr(data_tensor.DataTensor, 'charge_states'):
        title = '%s +%i Timepoint %s' % (
            data_tensor.DataTensor.name, data_tensor.DataTensor.charge_states[0], data_tensor.DataTensor.timepoint_idx)
    elif hasattr(data_tensor.DataTensor, 'charge_state'):
        title = '%s +%i Timepoint %s' % (
            data_tensor.DataTensor.name, data_tensor.DataTensor.charge_state, data_tensor.DataTensor.timepoint_idx)
    else:
        title = '%s Timepoint %s' % (data_tensor.DataTensor.name, data_tensor.DataTensor.timepoint_idx)

    plot_factor_data(retention_labels=data_tensor.DataTensor.retention_labels,
                     drift_labels=data_tensor.DataTensor.drift_labels,
                     mz_labels=data_tensor.DataTensor.mz_labels,
                     bins_per_isotope_peak=data_tensor.DataTensor.bins_per_isotope_peak,
                     tensor3=data_tensor.DataTensor.full_grid_out,
                     factors=data_tensor.DataTensor.factors,
                     output_path=output_path,
                     gauss_filter_params=data_tensor.gauss_params,
                     title=title,
                     from_dict=False)


def plot_factor_data_from_data_dict(factor_data, output_path):
    """

    :param factor_data:
    :param output_path:
    :return: None
    """

    title = ''
    if 'charge_states' in factor_data:
        title = '%s +%i Timepoint %s' % (
        factor_data['name'], factor_data['charge_states'][0], factor_data['timepoint_index'])
    elif 'charge_state' in factor_data:
        title = '%s +%i Timepoint %s' % (
            factor_data['name'], factor_data['charge_state'], factor_data['timepoint_index'])
    else:
        title = '%s Timepoint %s' % (factor_data['name'], factor_data['timepoint_index'])

    plot_factor_data(retention_labels=factor_data['retention_labels'],
                     drift_labels=factor_data['drift_labels'],
                     mz_labels=factor_data['mz_labels'],
                     bins_per_isotope_peak=factor_data['bins_per_isotope_peak'],
                     tensor3=factor_data['tensor_3d_grid'],
                     factors=factor_data['factors'],
                     output_path=output_path,
                     gauss_filter_params=factor_data['gauss_params'],
                     title=title,
                     from_dict=True)


def plot_factor_data_from_data_dict_file(factor_data_filepath, output_path=None):
    """
    plot factor data from factor data file .factor
    :param factor_data_filepath: .factor filepath
    :param output_path: output path
    :return: None. saves the figure
    """
    import _pickle as cpickle
    import zlib

    factor_data = cpickle.loads(zlib.decompress(open(factor_data_filepath, 'rb').read()))

    if output_path == None:
        output_path = factor_data_filepath + '.pdf'

    plot_factor_data_from_data_dict(factor_data=factor_data, output_path=output_path)


def plot_factor_from_commandline():
    """
    generate parser arguments for commandline execution
    :return:
    """
    import argparse

    parser = argparse.ArgumentParser(description='Plot factor data from a given .factor data file')
    parser.add_argument('--factor_data_path', action='store', help='file path to .factor file')
    parser.add_argument('--plot_output_path', action='store', help='output path for plot .pdf')

    args = parser.parse_args()

    plot_factor_data_from_data_dict_file(factor_data_filepath=args.factor_data_path,
                                         output_path=args.plot_output_path)

    return args


if __name__ == '__main__':

    plot_factor_from_commandline()

    # # to run it by modifying the script
    # factor_file_path = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/factor.factor'
    # output_path = '/Users/smd4193/Documents/MS_data/2021_lib15_ph6/factor.factor.pdf'
    # plot_factor_data_from_data_dict_file(factor_file_path, output_path=output_path)
