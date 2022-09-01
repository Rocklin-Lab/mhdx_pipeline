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
The first snakefile of the HDX_LIMIT pipeline. 
Determines signals to analyze by similarity to library protein expected values.
"""

import glob

def get_mem_mb(wildcards, attempt):
    return attempt * 2000

configfile: "config/config.yaml"

# Make flat list of all MS datafiles.
all_timepoint_files = []
for key in config["timepoints"]:
    for file in config[key]:
        all_timepoint_files.append(file)
        

if config['lockmass']:
    rule all:
        """
        Defines final outputs desired by pipeline run.
        """
        input:
            "resources/4_library_info/library_info.json",
            expand("resources/0_calibration/{mzml}_mz_calib_dict.pk", mzml=all_timepoint_files),
else:
    rule all:
            """
            Defines final outputs desired by pipeline run.
            """
            input:
                "resources/4_library_info/library_info.json",

if config['lockmass']:
    rule calibration_from_lockmass_0:
        input:
            "resources/2_mzml_gz/{mzml}.gz",
            "config/config.yaml"
        output:
            "resources/0_calibration/{mzml}_mz_calib_dict.pk",
            "results/plots/preprocessing/0_calibration/{mzml}_extracted_signals.pdf",
            "results/plots/preprocessing/0_calibration/{mzml}_degrees.pdf",
            "results/plots/preprocessing/0_calibration/{mzml}_kdes.pdf"
        benchmark:
            "results/benchmarks/0_calibration_from_lockmass.{mzml}.benchmark.txt"
        resources: mem_mb=get_mem_mb
        priority: 2
        script:
            "../scripts/hdx_limit/hdx_limit/preprocessing/0_calibration.py"

if config['lockmass']:
    rule read_imtbx_1:
        """
        Reads the identified peaks from the IMTBX .peaks.isotopes files made from undeuterated .mzML files. 
        Determines associations between identified peaks and expected masses of library proteins. 
        Feeds forward identified peaks with mass suffuciently similar to some library protein for further consideration. 
        """
        input:
            # The .peaks.isotopes files must be made in windows, but only for undeuterated MS runs.
            "config/config.yaml",
            "resources/0_isotopes/{undeut_fn}.peaks.isotopes",
            config["names_and_seqs"],
            "resources/0_calibration/{undeut_fn}_mz_calib_dict.pk",
        output:
            "resources/1_imtbx/{undeut_fn}_intermediate.csv",
            "results/plots/preprocessing/1_imtbx/{undeut_fn}_original_mz.pdf",
            "results/plots/preprocessing/1_imtbx/{undeut_fn}_adjusted_mz.pdf",
            "resources/1_imtbx/{undeut_fn}_mz_calib_dict.pk"
        benchmark:
            "results/benchmarks/1_read_imtbx.{undeut_fn}.benchmark.txt"
        resources: mem_mb=get_mem_mb
        script:
            "../scripts/hdx_limit/hdx_limit/preprocessing/1_imtbx_reader.py"
else:
    rule read_imtbx_1:
        """
        Reads the identified peaks from the IMTBX .peaks.isotopes files made from undeuterated .mzML files.
        Determines associations between identified peaks and expected masses of library proteins.
        Feeds forward identified peaks with mass suffuciently similar to some library protein for further consideration.
        """
        input:
            # The .peaks.isotopes files must be made in windows, but only for undeuterated MS runs.
            "config/config.yaml",
            "resources/0_isotopes/{undeut_fn}.peaks.isotopes",
            config["names_and_seqs"],
        output:
            "resources/1_imtbx/{undeut_fn}_intermediate.csv",
            "results/plots/preprocessing/{undeut_fn}_original_mz.pdf",
            "results/plots/preprocessing/{undeut_fn}_adjusted_mz.pdf",
            "resources/1_imtbx/{undeut_fn}_mz_calib_dict.pk"
        benchmark:
            "results/benchmarks/1_read_imtbx.{undeut_fn}.benchmark.txt"
        resources: mem_mb=get_mem_mb
        script:
            "../scripts/hdx_limit/hdx_limit/preprocessing/1_imtbx_reader.py"

rule gzip_mzmls_2:
    """
    Uses the pymzml module to make all .mzML files into randomly accessible .mzML.gz file, removes .mzMLs after conversion to save space.
    """
    input:
        "resources/0_mzml/{mzml}",
    output:
        "resources/2_mzml_gz/{mzml}.gz",
    resources: mem_mb=get_mem_mb
    shell:
        "python workflow/scripts/hdx_limit/hdx_limit/preprocessing/2_gzip_mzml.py {input} --delete_source --out_path {output}" 


rule make_ims_mz_tics_3:
    """
    Calculates total ionic current of an MS run at each LC retention timepoint, to be used by make_master_list.py
    """
    input:
        "resources/2_mzml_gz/{mzml}.gz",
    output:
        "resources/3_tics/{mzml}.ims.mz.tic.cpickle.zlib",
        "resources/3_tics/{mzml}_sum.txt"
    priority: 1
    resources: mem_mb=get_mem_mb
    shell:
        "python workflow/scripts/hdx_limit/hdx_limit/preprocessing/3_make_ims_mz_tics.py {input} --out_path {output[0]} --mzml_sum_outpath {output[1]}"


rule make_library_master_list_4:
    """
    Reads all candidate peaks from imtbx_reader output, 
    combines redundant references, clusters by retention-time,
    and generates final list of signals to consider for analysis.
    """
    input:
        config["names_and_seqs"],
        "config/config.yaml",
        expand("resources/2_mzml_gz/{mzml}.gz", mzml=all_timepoint_files[0]), # Pick a single mzML.gz
        expand(
            "resources/1_imtbx/{undeut_fn}_intermediate.csv", undeut_fn=config[0]
        ),
        glob.glob("resources/1_ms-unlimited/*_intermediate.csv"), # Get all intermediate files - imtbx and from unlimited
        expand("resources/3_tics/{mzml}.ims.mz.tic.cpickle.zlib", mzml=all_timepoint_files),
        expand("resources/3_tics/{mzml}_sum.txt", mzml=all_timepoint_files)
    output:
        "resources/4_library_info/library_info.json",
        "results/plots/preprocessing/4_make_library_master_list/stretched_times_plots.png",
        "resources/4_library_info/normalization_factors.csv",
        "results/plots/preprocessing/4_make_library_master_list/normalization_factors_plot.png",
        "results/plots/preprocessing/4_make_library_master_list/rt_correlation_plot.pdf"
    resources: mem_mb=get_mem_mb
    benchmark:
        "results/benchmarks/4_make_library_master_list.benchmark.txt"
    script:
        "../scripts/hdx_limit/hdx_limit/preprocessing/4_make_library_master_list.py"
