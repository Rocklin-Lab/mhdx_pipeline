"""
The third and final (maybe) snakefile of the HDX_LIMIT pipeline.
Extracts signals for quality-checked charge states, deconvolutes and outputs candidate IsotopeCluster objects.
Candidate IsotopeClusters from all charge states of an rt-group are pooled and a best-estimate HDX mass-addition timeseries is generated.
(Soon) Each mass-addition timeseries is used to estimate a set of amide hydrogen-deuterium exchange rates for the rt-group. 
(Soon) HDX rate sets are used to determine a Delta_G_unfolding for the protein. 
"""
import glob
import shutil
import pandas as pd
from collections import OrderedDict

def get_mem_mb(wildcards, attempt):
    return attempt * 2000

def optimize_paths_inputs(name, library_info):
    """Generate inputs to optimize_paths rule with rt-group name wildcard

        Args:
            name (str): The rt-group name passed as a wildcard.
            library_info (Pandas DataFrame): checked_library_info.json opened into a Pandas DF.

        Returns:
            name_inputs (list of strings): List of paths/to/input/files needed in optimize_paths.

    """
    name_inputs = []
    for key in config["timepoints"]:
        if len(config[key]) > 1:
            for charge in library_info.loc[library_info["name"]==name]['charge'].values:
                for file in config[key]:
                    name_inputs.append(
                    f"resources/9_subtensor_ics/{name}/{name}_charge{charge}_{file}.gz.cpickle.zlib"
                    )
        else:
            file = config[key][0]
            for charge in library_info.loc[library_info["name"]==name]['charge'].values:
                    name_inputs.append(
                    f"resources/9_subtensor_ics/{name}/{name}_charge{charge}_{file}.gz.cpickle.zlib"
                )
    return name_inputs


configfile: "config/config.yaml"
hdx_limit_dir = config["hdx_limit_dir"]

# Reads post idotp_check library_info.
library_info_fn = "resources/7_idotp_filter/checked_library_info.json"
library_info = pd.read_json(library_info_fn)

names = list(OrderedDict.fromkeys(library_info["name_rt-group"].values).keys()) # This is the Python-native version of an ordered set operation.
# Remove decoys from compute intensive factorization and path generation
if not config["decoys"]:
    names = [name for name in names if "decoy" not in name]
# Makes two zippable lists that are used for extract_tensors: repeated rt_group_names and their corresponding charges in order.
zippable_names = list(library_info["name_rt-group"].values)
zippable_charges = list(library_info["charge"].values)

# Creates three zippable lists that are used for mv_passing_tensors: rt-group names, charges, and undeut_mzmls.
mv_passing_tensors_zippable_names = []
mv_passing_tensors_zippable_charges = []
mv_passing_tensors_zippable_undeut_mzmls = []
for name, charge in zip(zippable_names, zippable_charges):
    for undeut_mzml in config[0]:
        mv_passing_tensors_zippable_names.append(name)
        mv_passing_tensors_zippable_charges.append(charge)
        mv_passing_tensors_zippable_undeut_mzmls.append(undeut_mzml)

mzmls = []
for tp in config["timepoints"]:
    for mzml in config[tp]:
        mzmls.append(mzml)


rule all:
    """
    Defines final outputs desired by pipeline run.
    """
    input:
        expand("resources/10_ic_time_series/{name}/monobody/{name}_winner_monobody.cpickle.zlib", name=names),
        expand("resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody.cpickle.zlib", name=names),
        expand("results/plots/ic_time_series/ajf_plots/multibody/{name}.pdf", name=names),
        expand("results/plots/ic_time_series/ajf_plots/monobody/{name}.pdf", name=names),
        "results/computational_resources_summary.pdf"


rule mv_passing_tensors_8:
    """
    Moves extracted undeuterated tensors that passed the idotp_check into a new working directory to limit redundancy and simplify input calling.
    """
    input:
        "config/config.yaml",
        library_info_fn,
        sorted(
            expand(
                "resources/5_tensors/{name}/{name}_charge{charge}_{mzml}.gz.cpickle.zlib",
                zip,
                name=mv_passing_tensors_zippable_names,
                charge=mv_passing_tensors_zippable_charges,
                mzml=mv_passing_tensors_zippable_undeut_mzmls
            )
        )
    output:
        sorted(
            expand(
                "resources/8_passing_tensors/{name}/{name}_charge{charge}_{mzml}.gz.cpickle.zlib",
                zip,
                name=mv_passing_tensors_zippable_names,
                charge=mv_passing_tensors_zippable_charges,
                mzml=mv_passing_tensors_zippable_undeut_mzmls
            )
        )
    resources: mem_mb=get_mem_mb
    benchmark:
        "results/benchmarks/8_mv_passing_tensors.benchmark.txt"
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/8_mv_passing_tensors.py"


rule extract_tensors_9:
    """
    Extract all identified tensors from each .mzML.gz.
    """
    input:
        library_info_fn,
        "resources/2_mzml_gz/{mzml}.gz",
        "config/config.yaml",
    output:
        expand(
            "resources/8_passing_tensors/{name}/{name}_charge{charge}_{{mzml}}.gz.cpickle.zlib",
            zip,
            name=zippable_names,
            charge=zippable_charges
        )
    resources: mem_mb=get_mem_mb
    params:
        use_rtdt_recenter=config["use_rtdt_recenter"]
    benchmark:
        "results/benchmarks/9_extract_tensors.{mzml}.gz.benchmark.txt"
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/5_extract_timepoint_tensors.py"


rule generate_tensor_ics_10:
    """
    Deconvolutes tensor signal and outputs candidate IsotopeClusters for path optimization.
    """
    input:
        "config/config.yaml",
        library_info_fn,
        "resources/4_library_info/normalization_factors.csv",
        expand(
            "resources/8_passing_tensors/{{name}}/{{name}}_charge{{charge}}_{mzml}.gz.cpickle.zlib",
            mzml=mzmls
        )
    output:
        expand(
            "resources/9_subtensor_ics/{{name}}/{{name}}_charge{{charge}}_{mzml}.gz.cpickle.zlib",
            mzml=mzmls
        )
    resources: mem_mb=get_mem_mb
    benchmark:
        "results/benchmarks/10_generate_tensor_ics.{name}/{name}_charge{charge}.benchmark.txt"
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/9_generate_tensor_ics.py"


rule generate_atcs_11:
    """
    Generate all timepoints cluster files
    """
    input:
        "config/config.yaml",
        lambda wildcards: optimize_paths_inputs(wildcards.name, library_info)
    output:
        "resources/10_ic_time_series/{name}/{name}_all_timepoint_clusters.cpickle.zlib"
    benchmark:
        "results/benchmarks/11_generate_atcs.{name}.benchmark.txt"
    resources: mem_mb=get_mem_mb
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/10_generate_atcs.py"


rule optimize_paths_tmp_12:
    """
    Takes all candidate ICs for all charges and timepoints of an rt-group and determines the best-estimate HDX mass-addition time series.
    """
    input:
        library_info_fn,
        "config/config.yaml",
        "resources/10_ic_time_series/{name}/{name}_all_timepoint_clusters.cpickle.zlib",
    output:
        "resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody_tmp.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_winner_scores_multibody_tmp.cpickle.zlib",
    params:
        rt_group_name = "{name}",
        tmp = True
    benchmark:
        "results/benchmarks/12_optimize_paths.{name}.benchmark.txt"
    resources: mem_mb=get_mem_mb
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/11_optimize_paths.py"


rule get_thresholds_13:
    input:
        expand("resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody_tmp.cpickle.zlib", name=names)
    output:
        "resources/10_ic_time_series/winners_tmp_dataframe.json",
        "resources/10_ic_time_series/thresholds.json"
    params:
        n_std = 2
    benchmark:
        "results/benchmarks/13_get_thresholds.benchmark.txt"
    resources: mem_mb=get_mem_mb
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/12_get_thresholds.py"


rule optimize_paths_14:
    """
    Takes all candidate ICs for all charges and timepoints of an rt-group and determines the best-estimate HDX mass-addition time series.
    """
    input:
        library_info_fn,
        "config/config.yaml",
        "resources/10_ic_time_series/{name}/{name}_all_timepoint_clusters.cpickle.zlib",
        "resources/10_ic_time_series/thresholds.json",
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
    params:
        rt_group_name = "{name}",
        tmp = False
    benchmark:
        "results/benchmarks/12_optimize_paths.{name}.benchmark.txt"
    resources: mem_mb=get_mem_mb
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/11_optimize_paths.py"


rule ajf_plot_14:
    input:
        "config/config.yaml",
        "resources/10_ic_time_series/{name}/{name}_all_timepoint_clusters.cpickle.zlib",
        "resources/10_ic_time_series/{name}/{name}_prefiltered_ics.cpickle.zlib",
        "resources/10_ic_time_series/{name}/multibody/{name}_winner_multibody.cpickle.zlib",
        "resources/10_ic_time_series/{name}/monobody/{name}_winner_monobody.cpickle.zlib"
    output:
        "results/plots/ic_time_series/ajf_plots/multibody/{name}.pdf",
        "results/plots/ic_time_series/ajf_plots/monobody/{name}.pdf"
    benchmark:
        "results/benchmarks/13_ajf_plots.{name}.benchmark.txt"
    resources: mem_mb=get_mem_mb
    shell:
         "python {hdx_limit_dir}/hdx_limit/pipeline/13_ajf_plot.py -c {input[0]} -a {input[1]} -f {input[2]} -w_multi {input[3]} -w_mono {input[4]} -o_multi {output[0]} -o_mono {output[1]}"

rule computational_resources:
    input:
        "results/benchmarks",
        expand("results/benchmarks/13_ajf_plots.{name}.benchmark.txt", name=names)
    output:
        "results/computational_resources_summary.pdf"
    resources: mem_mb=get_mem_mb
    shell:
        "python {hdx_limit_dir}/hdx_limit/core/generate_benchmark_plots.py -i {input[0]} -o {output[0]}"
