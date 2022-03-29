configfile: "config/config.yaml"
import glob
import shutil
import pandas as pd
from collections import OrderedDict

# Reads post idotp_check library_info.
library_info_fn = "resources/7_idotp_filter/checked_library_info.json"
library_info = pd.read_json(library_info_fn)

if config["use_rtdt_recenter"]:
    names = list(OrderedDict.fromkeys(library_info["name_recentered"].values).keys()) # This is the Python-native version of an ordered set operation.
    # Makes two zippable lists that are used for extract_tensors: repeated rt_group_names and their corresponding charges in order.
    zippable_names = list(library_info["name_recentered"].values)
    zippable_charges = list(library_info["charge"].values)
else:

    names = list(OrderedDict.fromkeys(library_info["name"].values).keys()) # This is the Python-native version of an ordered set operation.

    # Makes two zippable lists that are used for extract_tensors: repeated rt_group_names and their corresponding charges in order.
    zippable_names = list(library_info["name_recentered"].values)
    zippable_charges = list(library_info["charge"].values)


names = [name for name in names if name in [i.split('/')[-1] for i in glob.glob('resources/10_ic_time_series/*')] ]

rule all:
    """
    Defines final outputs desired by pipeline run.
    """
    input:
        expand("resources/10_ic_time_series/{name}/monobody/{name}_winner_monobody.cpickle.zlib", name=names),
        expand("resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody.cpickle.zlib", name=names),
        expand("results/plots/ic_time_series/ajf_plots/multibody/{name}.pdf", name=names),
        expand("results/plots/ic_time_series/ajf_plots/monobody/{name}.pdf", name=names)

rule optimize_paths_12:
    """
    Takes all candidate ICs for all charges and timepoints of an rt-group and determines the best-estimate HDX mass-addition time series.
    """
    input:
        library_info_fn,
        "config/config.yaml",
        "resources/10_ic_time_series/{name}/{name}_all_timepoint_clusters.cpickle.zlib"
    output:
        "resources/10_ic_time_series/{name}/{name}_prefiltered_ics.cpickle.zlib",
        "results/plots/ic_time_series/winner_plots/monobody/{name}_winner_path_monobody.pdf",
        "resources/10_ic_time_series/{name}/monobody/{name}_winner_monobody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/monobody/{name}_runners_monobody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/monobody/{name}_undeut_grounds_monobody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/monobody/{name}_winner_scores_monobody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/monobody/{name}_rtdt_com_cvs_monobody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/monobody/{name}_winner_monobody.cpickle.zlib.csv",
        "results/plots/ic_time_series/winner_plots/multibody/{name}_winner_path_multibody.pdf",
        "resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_runners_multibody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_undeut_grounds_multibody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_winner_scores_multibody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_rtdt_com_cvs_multibody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody.cpickle.zlib.csv",
#         "results/plots/ic_time_series/ajf_plot/{name}.pdf"
    params:
        rt_group_name = "{name}"
    benchmark:
        "results/benchmarks/12_optimize_paths.{name}.benchmark.txt"
    conda:
        "../envs/full_hdx_env.yml"
    script:
        "../scripts/hdx_limit/hdx_limit/pipeline/11_optimize_paths.py"

rule ajf_plot_13:
    input:
        "config/config.yaml",
        "resources/10_ic_time_series/{name}/{name}_all_timepoint_clusters.cpickle.zlib",
        "resources/10_ic_time_series/{name}/{name}_prefiltered_ics.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/monobody/{name}_winner_monobody.cpickle.zlib"
    output:
        "results/plots/ic_time_series/ajf_plots/multibody/{name}.pdf",
        "results/plots/ic_time_series/ajf_plots/monobody/{name}.pdf"
    script:
         "../scripts/hdx_limit/hdx_limit/pipeline/12_ajf_plot.py"
# TEST
