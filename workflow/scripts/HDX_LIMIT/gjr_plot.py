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
import zlib
import numpy as np
import _pickle as cpickle
from matplotlib import pyplot as plt


def plot_gjr_(winner, undeut_grounds, output_path, prefix='winner_plot'):
    """
    plot path output given the winner ic list and undeut grounds list
    :param winner: winner list
    :param undeut_grounds: undeut grounds list
    :param output_path: plot output path .pdf
    :param prefix: used in plot title
    :return:
    """
    fig = plt.figure(figsize=(20, 24))

    protname = prefix  # 'EHEE_rd1_0284.pdb_5.73355'
    idotp = undeut_grounds[1][winner[0].charge_states[0]]
    for i, x in enumerate(winner):

        ##########
        #
        # INTEGRATED MZ
        #
        ##########

        ax = plt.subplot(len(winner), 3, (3 * i) + 1)
        plt.plot(x.baseline_integrated_mz / max(x.baseline_integrated_mz), linewidth=1)
        plt.scatter(range(len(x.baseline_integrated_mz)), x.baseline_integrated_mz / max(x.baseline_integrated_mz), s=5)
        plt.yticks([])
        plt.xticks(range(0, len(x.baseline_integrated_mz) + 5, 5))
        ax.set_xticklabels(range(0, len(x.baseline_integrated_mz) + 5, 5), fontsize=8)
        ax.tick_params(length=3, pad=3)
        plt.grid(axis='x', alpha=0.25)
        plt.ylim(0, 1.4)
        plt.plot([x.baseline_integrated_mz_com, x.baseline_integrated_mz_com], [1.0, 1.1], color='black')

        plt.text(
            x.baseline_integrated_mz_com,
            1.1,
            "%.1f" % (x.baseline_integrated_mz_com),
            horizontalalignment="center",
            verticalalignment="bottom",
            fontsize=8
        )

        # for runner in runners[i]:
        #    if (
        #        runner.dt_ground_fit > x.dt_ground_fit
        #        and runner.rt_ground_fit > x.rt_ground_fit
        #        and abs(
        #            runner.log_baseline_auc
        #            - undeut_grounds[0][runner.charge_states[0]].log_baseline_auc
        #        )
        #        < abs(
        #            x.log_baseline_auc
        #            - undeut_grounds[0][x.charge_states[0]].log_baseline_auc
        #        )
        #    ):
        #
        #        plt.plot(
        #            runner.baseline_integrated_mz / max(runner.baseline_integrated_mz)
        #        )

        idotp_q = dict((v, k) for k, v in undeut_grounds[1].items())
        max_idotp = max(idotp_q.keys())

        if i == 0:
            plt.title("%s idotp %.3f" % (protname, idotp))
            plt.text(
                1.0,
                0.1,
                "best idotp +%.0f %.3f" % (idotp_q[max_idotp], max_idotp),
                horizontalalignment="right",
                verticalalignment="bottom",
                transform=ax.transAxes,
            )

        add_delta_com_text = ''
        if i > 0:
            delta_mz_com = (x.baseline_integrated_mz_com - winner[i - 1].baseline_integrated_mz_com)
            add_delta_com_text = '\nd_mz %.1f' % (delta_mz_com)

        plt.text(
            1.0,
            1.0,
            "+%.0f, %s factors%s" % (x.charge_states[0], x.n_factors, add_delta_com_text),
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
        )

        plt.text(
            0.01,
            1.2,
            '%s t %s' % (i, x.timepoint_idx),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ##########
        #
        # RAW ISOTOPIC DISTRIBUTION
        #
        ##########

        reshape_mz = np.reshape(x.baseline_subtracted_mz, (-1, x.bins_per_isotope_peak))
        reshape_mz_labels = np.reshape(x.mz_labels, (-1, x.bins_per_isotope_peak))
        nonzero_isotopes = np.nonzero(np.sum(reshape_mz, axis=1))[0]

        xtick_lims = [reshape_mz_labels[nonzero_isotopes[0]][int(x.bins_per_isotope_peak / 2)],
                      reshape_mz_labels[nonzero_isotopes[-1]][int(x.bins_per_isotope_peak / 2)]]

        plot_raw_mz = []
        plot_raw_mz_labels = []
        mz_label_gap = x.mz_labels[1] - x.mz_labels[0]
        for j in nonzero_isotopes:
            plot_raw_mz += ([0] + list(reshape_mz[j]) + [0])
            plot_raw_mz_labels += ([reshape_mz_labels[j][0] - mz_label_gap] + list(reshape_mz_labels[j]) + [
                reshape_mz_labels[j][-1] + mz_label_gap])

        ax = plt.subplot(len(winner), 6, (6 * i) + 3)
        plt.plot(plot_raw_mz_labels, plot_raw_mz / max(plot_raw_mz))
        plt.yticks([])
        plt.xticks(np.linspace(xtick_lims[0], xtick_lims[1], 3))
        ax.set_xticklabels(['%.4f' % x for x in np.linspace(xtick_lims[0], xtick_lims[1], 3)], fontsize=8)
        if i == 0:
            plt.title("Raw isotopic distribution")
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.text(
            0.05,
            1.0,
            "gauss_rms %.2f nearest %.2f" % (x.baseline_integrated_mz_rmse, x.nearest_neighbor_correlation),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )
        plt.ylim(0, 1.4)

        ###########
        #
        # PEAK ERROR
        #
        ###########

        ax = plt.subplot(len(winner), 6, (6 * i) + 4)
        for j in nonzero_isotopes:
            new_x = np.linspace(0, 1, x.bins_per_isotope_peak)
            plt.plot(new_x, reshape_mz[j])
            plt.scatter(new_x, reshape_mz[j], s=20)
        ylim = ax.get_ylim()
        plt.plot([0.5, 0.5], ylim, color="black")
        plt.ylim([ylim[0], ylim[1] * 1.25])

        plt.yticks([])
        plt.xticks([])
        if i == 0:
            plt.title("Overlapped peaks")
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # peak error text
        plt.text(
            0.05,
            1.0,
            "centering err %.2f" % (x.peak_error),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )

        ##########
        #
        # RETENTION
        #
        ##########

        ax = plt.subplot(len(winner), 9, (9 * i) + 7)
        plt.plot(x.rts / max(x.rts))
        plt.plot(
            undeut_grounds[0][x.charge_states[0]].rts
            / max(undeut_grounds[0][x.charge_states[0]].rts),
            color="red",
        )
        plt.yticks([])
        plt.xticks([0, len(x.rts) / 2, len(x.rts)])
        ax.tick_params(length=3, pad=3)
        ax.set_xticklabels(['%.2f' % x for x in np.linspace(x.retention_labels[0], x.retention_labels[-1], 3)],
                           fontsize=8)

        xTick_objects = ax.xaxis.get_major_ticks()
        xTick_objects[0].label1.set_horizontalalignment('left')  # left align first tick
        xTick_objects[-1].label1.set_horizontalalignment('right')  # right align last tick

        # xlabels = [tick.label for tick in ax.xaxis.get_major_ticks()]

        # for i in range(3):
        #    xlabels[i].set_horizontalalignment(alignments[i])
        #
        #    pass

        plt.text(
            0.05,
            1.0,
            "err %.1f fit %.2f" % (x.rt_ground_err, x.rt_ground_fit),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )
        plt.ylim(0, 1.4)
        if i == 0:
            plt.title("RT (red=undeut)")

        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ##########
        #
        # DRIFT
        #
        ##########

        ax = plt.subplot(len(winner), 9, (9 * i) + 8)
        plt.plot(x.dts / max(x.dts))
        plt.plot(
            undeut_grounds[0][x.charge_states[0]].dts
            / max(undeut_grounds[0][x.charge_states[0]].dts),
            color="red",
        )

        plt.xticks([0, len(x.dts) / 2, len(x.dts)])
        ax.tick_params(length=3, pad=3)
        ax.set_xticklabels(['%.2f' % x for x in np.linspace(x.drift_labels[0], x.drift_labels[-1], 3)], fontsize=8)

        xTick_objects = ax.xaxis.get_major_ticks()
        xTick_objects[0].label1.set_horizontalalignment('left')  # left align first tick
        xTick_objects[-1].label1.set_horizontalalignment('right')  # right align last tick

        plt.text(
            0.05,
            1.0,
            "err %.1f fit %.2f" % (x.dt_ground_err, x.dt_ground_fit),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )
        plt.yticks([])

        plt.ylim(0, 1.4)
        if i == 0:
            plt.title("Drift time")
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ##########
        #
        # TOTAL INTENSITY
        #
        ##########

        ax = fig.add_subplot(len(winner), 9, (9 * i) + 9)
        plt.plot(x.dts / max(x.dts))
        plt.bar(
            [0, 1],
            [
                sum(undeut_grounds[0][x.charge_states[0]].baseline_subtracted_mz)
                * undeut_grounds[0][x.charge_states[0]].outer_rtdt,
                sum(x.baseline_subtracted_mz) * x.outer_rtdt,
            ],
            color=["red", "blue"],
        )
        plt.yticks([])
        plt.xticks([])
        plt.xlim(-0.5, 5)
        plt.text(
            1.0,
            1.0,
            "%.1e"
            % (
                sum(undeut_grounds[0][x.charge_states[0]].baseline_subtracted_mz)
                * undeut_grounds[0][x.charge_states[0]].outer_rtdt
            ),
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
            color="red",
        )
        plt.text(
            1.0,
            1.0,
            "\n%.1e" % (sum(x.baseline_subtracted_mz) * x.outer_rtdt),
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
            color="blue",
        )

        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if i == 0:
            plt.title("Magnitude")

    plt.subplots_adjust(hspace=0.5, wspace=0.1, top=0.95)
    plt.savefig(output_path, bbox_inches="tight")


def plot_gjr_from_files(winner_path_fpath, undeut_grounds_fpath, output_path, prefix='winner_plot'):
    winner = cpickle.loads(zlib.decompress(open(winner_path_fpath, "rb").read()))

    undeut_grounds = cpickle.loads(zlib.decompress(open(undeut_grounds_fpath, "rb").read()))

    plot_gjr_(winner=winner, undeut_grounds=undeut_grounds, output_path=output_path, prefix=prefix)


def plot_path_from_commandline():
    """
    generate parser arguments for commandline execution
    :return:
    """
    import argparse

    parser = argparse.ArgumentParser(description='Plot factor data from a given .factor data file')
    parser.add_argument('--winners_path', '-w', action='store', help='file path to winner.cpickle.zlib')
    parser.add_argument('--undeut_grounds_path', '-u', action='store', help='file path to undeut_grounds.cpickle.zlib')
    parser.add_argument('--prefix', '-p', action='store', help='prefix str to add to the plot', default='winner_plot')
    parser.add_argument('--plot_output_path', '-o', action='store', help='output path for plot .pdf')

    args = parser.parse_args()

    plot_gjr_from_files(winner_path_fpath=args.winners_path,
                        undeut_grounds_fpath=args.undeut_grounds_path,
                        output_path=args.prefix,
                        prefix=args.plot_output_path)

    return args


if __name__ == '__main__':

    plot_path_from_commandline()
