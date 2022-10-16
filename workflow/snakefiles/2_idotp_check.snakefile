import glob
import pandas as pd
from collections import OrderedDict

def get_mem_mb(wildcards, attempt):
    return attempt * 2000

def idotp_check_inputs(config, rt_group_name, charge):
    """Wildcard-based input-file-path generator for idotp_check rule, writes based on rt_group and charge.

    Args:
        rt_group_name (str): Value in 'name' field of library_info.json, shared by groups of mass-agreeing signals close in rt.
        charge (int): Net charge of signal, unique for each member of an rt-group.

    Returns:
        inputs (list of strs): List of paths to extracted tensors from undeuterated timepoints for a given rt-group charge state.

    """
    inputs = []
    if len(config[0]) > 1:
        for undeut_fn in config[0]:
            inputs.append(
                "resources/5_tensors/" + rt_group_name + "/" + rt_group_name + "_" + "charge" + str(charge) + "_" + undeut_fn + ".gz.cpickle.zlib"
            )
    else: # Handle 0th timepoint with no replicates.
        undeut_fn = config[0][0]
        inputs.append(
            "resources/5_tensors/" + rt_group_name + "/" + rt_group_name + "_" + "charge" + str(charge) + "_" + undeut_fn + ".gz.cpickle.zlib"
        )
    return inputs

configfile: "config/config.yaml"

hdx_limit_dir = config["hdx_limit_dir"]

# Read list of candidate POI charge states produced by preprocessing snakefile
library_info_fn = "resources/4_library_info/library_info.json"
library_info = pd.read_json(library_info_fn)
names = list(OrderedDict.fromkeys(library_info["name"].values).keys()) # This is the Python-native version of an ordered set operation.

# Makes two zippable lists: repeated rt_group_names and their corresponding charges in order, 
# used for extract_tensors and idotp_filter rules.
zippable_names = list(library_info["name"].values)
zippable_charges = list(library_info["charge"].values)

# Make flat list of replicate .mzML files in timepoint order.
mzml_list = []
for timepoint in config["timepoints"]:
    for fn in config[timepoint]:
        mzml_list.append(fn)

rule all:
    """
    Defines final outputs desired by pipeline run.
    """
    input:
        "resources/7_idotp_filter/checked_library_info.json"


rule extract_tensors_5:
    """
    Extract all identified tensors from each .mzML.gz.
    """
    input:
        library_info_fn,
        "resources/2_mzml_gz/{mzml}.gz",
        "config/config.yaml",
        "resources/0_calibration/{mzml}_mz_calib_dict.pk",
        "resources/1_imtbx/{mzml}_mz_calib_dict.pk"
    output:
        expand(
            "resources/5_tensors/{name}/{name}_charge{charge}_{{mzml}}.gz.cpickle.zlib",
            zip,
            name=zippable_names,
            charge=zippable_charges
        )
    params:
        use_rtdt_recenter=False
    resources: mem_mb=get_mem_mb
    benchmark:
        "results/benchmarks/5_extract_tensors.{mzml}.gz.benchmark.txt"
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/5_extract_timepoint_tensors.py"


rule idotp_check_6:
    """
    Test processed signal from each extracted tensor from undeuterated timepoints against the theoretical signal.
    """
    input:
        library_info_fn,
        "config/config.yaml",
        "resources/4_library_info/normalization_factors.csv",
        lambda wildcards: idotp_check_inputs(config,wildcards.name, wildcards.charge),
    output:
        "resources/6_idotp_check/{name}/{name}_charge{charge}_idotp_check.json",
        "resources/6_idotp_check/{name}/{name}_charge{charge}.cpickle.zlib",
    resources: mem_mb=get_mem_mb
    benchmark:
        "results/benchmarks/6_idotp_check.{name}_charge{charge}.benchmark.txt"
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/6_idotp_check.py"


rule idotp_filter_7:
    """
    Read results of idotp_check for all undeuterated tensors and apply a high-pass filter. 
    Produces a list of passing indices and and edited version of library_info. 
    """
    input:
        "config/config.yaml",
        library_info_fn,
        expand(
            "resources/6_idotp_check/{name}/{name}_charge{charge}_idotp_check.json",
            zip,
            name=zippable_names,
            charge=zippable_charges
        ),
        expand(
            "resources/6_idotp_check/{name}/{name}_charge{charge}.cpickle.zlib",
            zip,
            name=zippable_names,
            charge=zippable_charges
        )
    output:
        "resources/7_idotp_filter/checked_library_info.json",
        "results/plots/idotp_filter/idotp_distribution.png",
        "results/plots/idotp_filter/undeuderated_deviations.pdf",
        "results/plots/idotp_filter/fdr_stats.pdf"
    resources: mem_mb=get_mem_mb
    benchmark:
        "results/benchmarks/7_idotp_filter.benchmark.txt"
    script:
        f"{hdx_limit_dir}/hdx_limit/pipeline/7_idotp_filter.py"


