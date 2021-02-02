import os
import glob
import zlib
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import _pickle as cpickle
from matplotlib import pyplot as plt

sys.path.append(os.getcwd() + "/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx


def snakemake_makeplot(prefix, inputs, output):
    with open(inputs[0], "rb") as file:
        winner = cpickle.loads(zlib.decompress(file.read()))
    with open(inputs[1], "rb") as file:
        runners = cpickle.loads(zlib.decompress(file.read()))
    with open(inputs[2], "rb") as file:
        undeut_grounds = cpickle.loads(zlib.decompress(file.read()))

    plt.figure(figsize=(20, 22))

    protname = prefix  #'EHEE_rd1_0284.pdb_5.73355'
    idotp = undeut_grounds[1][winner[0].charge_states[0]]
    for i, x in enumerate(winner):

        ax = plt.subplot(len(winner), 3, (3 * i) + 1)
        plt.plot(x.baseline_integrated_mz / max(x.baseline_integrated_mz))
        plt.yticks([])
        plt.xticks([])
        sns.despine()

        for runner in runners[i]:
            if (
                runner.dt_ground_fit > x.dt_ground_fit
                and runner.rt_ground_fit > x.rt_ground_fit
                and abs(
                    runner.log_baseline_auc
                    - undeut_grounds[0][runner.charge_states[0]].log_baseline_auc
                )
                < abs(
                    x.log_baseline_auc
                    - undeut_grounds[0][x.charge_states[0]].log_baseline_auc
                )
            ):

                plt.plot(
                    runner.baseline_integrated_mz / max(runner.baseline_integrated_mz)
                )

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

        plt.text(
            1.0,
            1.0,
            "+%.0f, %s factors" % (x.charge_states[0], x.n_factors),
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
        )

        plt.text(
            0.01,
            1.0,
            str(i),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )

        # ax = plt.subplot(len(winner), 6, (6*i)+3)
        # plt.plot(x.factor_mz_data)
        # plt.yticks([])
        # plt.xticks([])
        # plt.text(1.0,1.0,'%s factors' % x.n_factors,
        #         horizontalalignment='right',verticalalignment='top',transform=ax.transAxes)
        # sns.despine()

        ax = plt.subplot(len(winner), 6, (6 * i) + 3)
        plt.plot(np.trim_zeros(x.baseline_subtracted_mz))
        plt.yticks([])
        plt.xticks([])
        if i == 0:
            plt.title("Raw isotopic distribution")
        sns.despine()

        ax = plt.subplot(len(winner), 6, (6 * i) + 4)
        for low, hi in zip(x.lows, x.highs):
            if sum(x.baseline_subtracted_mz[low:hi]) > 0:
                new_x = np.linspace(0, 1, hi + 1 - low)
                plt.plot(new_x, x.baseline_subtracted_mz[low : hi + 1])
                plt.scatter(new_x, x.baseline_subtracted_mz[low : hi + 1], s=20)
        ylim = ax.get_ylim()
        plt.plot([0.5, 0.5], ylim, color="black")
        plt.ylim(ylim)

        plt.yticks([])
        plt.xticks([])
        if i == 0:
            plt.title("Overlapped peaks")
        sns.despine()

        ax = plt.subplot(len(winner), 9, (9 * i) + 7)
        plt.plot(x.rts / max(x.rts))
        plt.plot(
            undeut_grounds[0][x.charge_states[0]].rts
            / max(undeut_grounds[0][x.charge_states[0]].rts),
            color="red",
        )
        plt.yticks([])
        plt.xticks([])
        plt.text(
            1.0,
            1.0,
            "%.3f" % x.rt_ground_fit,
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
        )
        plt.ylim(0, 1.4)
        if i == 0:
            plt.title("RT (red=undeut)")
        # sns.despine()

        ax = plt.subplot(len(winner), 9, (9 * i) + 8)
        plt.plot(x.dts / max(x.dts))
        plt.plot(
            undeut_grounds[0][x.charge_states[0]].dts
            / max(undeut_grounds[0][x.charge_states[0]].dts),
            color="red",
        )
        plt.text(
            1.0,
            1.0,
            "%.3f" % x.dt_ground_fit,
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
        )
        plt.yticks([])
        plt.xticks([])
        plt.ylim(0, 1.4)
        if i == 0:
            plt.title("Drift time")
        # sns.despine()

        ax = plt.subplot(len(winner), 9, (9 * i) + 9)
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
        sns.despine()
        if i == 0:
            plt.title("Magnitude")
    plt.savefig(output[0], bbox_inches="tight")


prefix = snakemake.input[0][3:].replace("_winner.cpickle.zlib", "")
snakemake_makeplot(snakemake.wildcards["name"], snakemake.input, snakemake.output)
