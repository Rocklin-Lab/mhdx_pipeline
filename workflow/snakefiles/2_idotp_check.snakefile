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

"""
The second Snakefile of the HDX_LIMIT-Pipeline. 
Checks quality of identified signals in undeuterated data 
and removes signals below user-defined quality thresholds from consideration. 
"""

configfile: "config/config.yaml" # Sets 'config' global object.
import glob
import pandas as pd
from collections import OrderedDict

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

def idotp_check_inputs(rt_group_name, charge):
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
    

rule all:
    """
    Defines final outputs desired by pipeline run.
    """
    input:
        "resources/7_idotp_filter/filter_passing_indices.csv",
        "resources/7_idotp_filter/checked_library_info.json"


if config['polyfit_calibration']:
    rule extract_tensors_5:
        """
        Extract all identified tensors from each .mzML.gz. 
        """
        input:
            library_info_fn,
            "resources/2_mzml_gz/{mzml}.gz",
            "config/config.yaml",
            expand(
                "results/1_imtbx/{undeut_fn}_mz_calib_dict.pk", undeut_fn=config[0][0]
            )
        output:
            expand(
                "resources/5_tensors/{name}/{name}_charge{charge}_{{mzml}}.gz.cpickle.zlib",
                zip,
                name=zippable_names,
                charge=zippable_charges
            )
        conda: 
            "envs/full_hdx_env.yml"
        benchmark:
            "results/benchmarks/5_extract_tensors.{mzml}.gz.benchmark.txt"

        script:
            "scripts/5_extract_timepoint_tensors.py"
else:
    rule extract_tensors_5:
        """
        Extract all identified tensors from each .mzML.gz. 
        """
        input:
            library_info_fn,
            "resources/2_mzml_gz/{mzml}.gz",
            "config/config.yaml"
        output:
            expand(
                "resources/5_tensors/{name}/{name}_charge{charge}_{{mzml}}.gz.cpickle.zlib",
                zip,
                name=zippable_names,
                charge=zippable_charges
            )
        conda: 
            "envs/full_hdx_env.yml"
        benchmark:
            "results/benchmarks/5_extract_tensors.{mzml}.gz.benchmark.txt"
        script:
            "scripts/5_extract_timepoint_tensors.py"


rule idotp_check_6:
    """
    Test processed signal from each extracted tensor from undeuterated timepoints against the theoretical signal. 
    """
    input:
        library_info_fn,
        "config/config.yaml",
        "resources/4_library_info/normalization_factors.csv",
        lambda wildcards: idotp_check_inputs(wildcards.name, wildcards.charge),
    output:
        "resources/6_idotp_check/{name}/{name}_charge{charge}_idotp_check.json",
        expand("resources/6_idotp_check/{{name}}/{{name}}_charge{{charge}}_{file}.cpickle.zlib.factor", file=config[0]),
        expand("results/plots/factors/{{name}}/{{name}}_charge{{charge}}_{file}.cpickle.zlib.factor.pdf", file=config[0])
    conda: 
        "envs/full_hdx_env.yml"
    benchmark:
        "results/benchmarks/6_idotp_check.{name}_charge{charge}.benchmark.txt"
    script:
        "scripts/6_idotp_check.py"


rule idotp_filter_7:
    """
    Read results of idotp_check for all undeuterated tensors and apply a high-pass filter. 
    Produces a list of passing indices and and edited version of library_info. 
    """
    input: 
        library_info_fn,
        expand(
            "resources/6_idotp_check/{name}/{name}_charge{charge}_idotp_check.json",
            zip,
            name=zippable_names,
            charge=zippable_charges
        )
    output:
        "resources/7_idotp_filter/filter_passing_indices.csv",
        "resources/7_idotp_filter/checked_library_info.json",
        "results/plots/idotp_distribution.png"
    conda: 
        "envs/full_hdx_env.yml"
    benchmark:
        "results/benchmarks/7_idotp_filter.benchmark.txt"
    script:
        "scripts/7_idotp_filter.py"

