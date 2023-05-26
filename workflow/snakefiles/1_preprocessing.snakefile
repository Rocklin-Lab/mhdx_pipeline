import os
import glob
from pathlib import Path
from hdx_limit.core.generate_decoy_database import main as generate_decoy_database

def get_mem_mb(wildcards, attempt):
    return attempt * 2000

configfile : "config/config.yaml"  # Sets 'config' global object.

hdx_limit_dir = config["hdx_limit_dir"]

# Generate decoy dataset
if not os.path.exists(f'{os.path.dirname(config["names_and_seqs"])}/decoys.csv'):
    generate_decoy_database(name_mass_seq=config["names_and_seqs"],
                            decoy_size=config["decoy_level"],
                            output_path=f'{os.path.dirname(config["names_and_seqs"])}/decoys.csv',
                            )

config["names_and_seqs"] = f'{os.path.dirname(config["names_and_seqs"])}/decoys.csv'


# Make flat list of all MS datafiles.
all_timepoint_files = []
for key in config["timepoints"]:
    for file in config[key]:
        all_timepoint_files.append(file)


rule all:
    input:
        "resources/4_library_info/library_info.json",
        expand("resources/0_calibration/{mzml}_mz_calib_dict.pk", mzml=all_timepoint_files) if config["lockmass"] else []

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
        f"{hdx_limit_dir}/hdx_limit/preprocessing/0_calibration.py"

rule read_imtbx_1:
    """
    Reads the identified peaks from the IMTBX .peaks.isotopes files made from undeuterated .mzML files.
    Determines associations between identified peaks and expected masses of library proteins.
    Feeds forward identified peaks with mass suffuciently similar to some library protein for further consideration.
    """
    input:
        "config/config.yaml",
        "resources/0_isotopes/{undeut_fn}.peaks.isotopes",
        config["names_and_seqs"],
        "resources/0_calibration/{undeut_fn}_mz_calib_dict.pk" if config["lockmass"] else []
    output:
        "resources/1_imtbx/{undeut_fn}_intermediate.csv",
        "results/plots/preprocessing/1_imtbx/{undeut_fn}_errors_kde.pdf",
    benchmark:
        "results/benchmarks/1_read_imtbx.{undeut_fn}.benchmark.txt"
    resources: mem_mb=get_mem_mb
    script:
         f"{hdx_limit_dir}/hdx_limit/preprocessing/1_imtbx_reader.py"

rule gzip_mzmls_2:
    """
    Uses the pymzml module to make all .mzML files into randomly accessible .mzML.gz file, removes .mzMLs after conversion to save space.
    """
    input:
        "resources/0_mzml/{mzml}",
    output:
        "resources/2_mzml_gz/{mzml}.gz",
    benchmark:
        "results/benchmarks/2_gzip_mzml.{mzml}.benchmark.txt"
    resources: mem_mb=get_mem_mb
    shell:
        "python {hdx_limit_dir}/hdx_limit/preprocessing/2_gzip_mzml.py {input} --delete_source --out_path {output}"

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
    benchmark:
        "results/benchmarks/3_make_ims_mz_tics.{mzml}.benchmark.txt"
    resources: mem_mb=get_mem_mb
    shell:
        "python {hdx_limit_dir}/hdx_limit/preprocessing/3_make_ims_mz_tics.py {input} --out_path {output[0]} --mzml_sum_outpath {output[1]}"

rule make_library_master_list_4:
    """
    Reads all candidate peaks from imtbx_reader output,
    combines redundant references, clusters by retention-time,
    and generates final list of signals to consider for analysis.
    """
    input:
        "config/config.yaml",
        glob.glob("resources/1_*/*_intermediate.csv"), # Get all intermediate files - imtbx and from unlimited
#        expand("resources/2_mzml_gz/{mzml}.gz", mzml=all_timepoint_files[0]), # Pick a single mzML.gz
#        expand(
#            "resources/1_imtbx/{undeut_fn}_intermediate.csv", undeut_fn=config[0]
#        ),
#        expand("resources/3_tics/{mzml}.ims.mz.tic.cpickle.zlib", mzml=all_timepoint_files),
        expand("resources/3_tics/{mzml}_sum.txt", mzml=all_timepoint_files)
    output:
        "resources/4_library_info/library_info.json",
        "resources/4_library_info/normalization_factors.csv",
        "results/plots/preprocessing/4_make_library_master_list/normalization_factors_plot.png",
        "results/plots/preprocessing/4_make_library_master_list/rt_correlation_plot.pdf",
        "results/plots/preprocessing/4_make_library_master_list/rt_distribution_plot.pdf",
        #"results/plots/preprocessing/4_make_library_master_list/stretched_times_plots.png",
    resources: mem_mb=get_mem_mb
    benchmark:
        "results/benchmarks/4_make_library_master_list.benchmark.txt"
    script:
        f"{hdx_limit_dir}/hdx_limit/preprocessing/4_make_library_master_list.py"
