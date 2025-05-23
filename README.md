# The mHDX-MS Pipeline with Snakemake

### Version 1.25.35

## Overview

The `mHDX-MS` pipeline is designed for processing mass spectrometry data from hydrogen-deuterium exchange (HDX-MS) experiments from LC-IMS-MS analysis. Utilizing `Snakemake`, it ensures scalable, reproducible data analysis. 

To facilitate reproducibility and user onboarding, we provide a set of curated example files associated with one of our experimental libraries. These files, available under the `Demo` section, include all necessary inputs to execute a full pipeline run, including protein sequences, processed mass spectrometry data, and a template configuration file for both pH6 and pH9 experiments. 

---


## Installing mhdx_tools

Step 1: Install Miniconda (if not already installed). Follow the instruction for your operating system: 

```
https://docs.anaconda.com/miniconda/install/#quick-command-line-install
```

Step 2: Download `mhdx_tools` package
```
git clone https://github.com/Rocklin-Lab/mhdx_tools.git
```

Step 3: Install `mamba` and create a new environment
```
conda install -n base -c conda-forge mamba
mamba env create -f mhdx_tools/environment.yaml
```

Step 4: Activate the environment and install `mhdx_tools`

```
conda activate mhdx
python -m pip install ./mhdx_tools
```


## Setting up a new experiment

For each mHDX-MS experiment, clone a new pipeline repository. For example, within `{library}/{date-of-experiment}_{pH6/pH9}/`

`git clone https://github.com/Rocklin-Lab/mhdx_pipeline.git`

The pipeline expects a set of input files:

1) a csv file contaning protein name, protein sequence, and monoisotopic mass. Copy this file to `resources/0_names_seqs_masses` <br />
e.g. names_and_seqs.csv <br />
```
name,sequence,mono_mass
PDB2PJV,HMAVGIGALFLGFLGAAGSTVGAASGGGKKKKK,3085.722267
PDB2N92,HMSWLSKTAKKLENSAKKRISEGIAIAIQGGPR,3604.9987800000013
PDB2LZP,HMDTEIIGGLTIPPVVALVVMSRFGFFAHLLPR,3632.972752
```

2) mzML files from HX-MS experiment. Current pipeline module to read .gz files was fully tested with `ms-convert` from `ProteomeWizard 3.0.21193-ccb3e0136`. <br />
Copy `mzML` files to `resources/0_mzml`. Preferably, the user can provide already `gziped` versions of the `mzML` files. In this case, copy `mzML.gz` files to `resources/2_mzml_gz`. 

Some options I suggest when converting files:
```bash
- Options
Binary encoding precision: 32-bit
Use zlib compression: True
Package in gzip: True
All other Options not selected

- Filters
Zero Samples: Remove --> Add
```

3) isotope file for each undeuterated sample. This file is generated with `IMTBX` and `Grppr`. Example commands are provided below. More information can be found at `https://dmtavt.com/IMTBX`. Copy these files to `resources/0_isotopes`. They should be named `{undeuterated-file-name}.mzML.peaks.isotopes`

### Running IMTBX and Grppr
The following commands run the `IMTBX` and `Grppr` executables to process peaks:

```bash
.\IMTBX.exe peaks \
    -vw \
    --mode ScanByScan \
    --clean True \
    --cut 10 \
    --filter 1 1 1 1 1 1 \
    --hyper 0.4 0.5 0.5 \
    -n True \
    --noise True \
    --noiseWnd 2 4 \
    --snr 0.5 \
    -i $file \
    --lock False \
    -o "+../Results/" \
    --ignore-warnings False \
    --orig

java -jar .\grppr-0.3.21.jar \
    -d 1 4 \
    --dSnr 1.0 1.0 \
    --isoAlg AVERAGINE \
    --isoMzV 30.0 \
    --isoMzT PPM \
    --isoTi 1.0 \
    -g False \
    --zLo 3 \
    --zHi 15 \
    -w \
    -v \
    -i $file
```

5) edit config/config.yaml to contain the specific paths to those files along with the set of timepoints for your experiment.

Import blocks to edit
```
"mhdx_tools_dir": <path-to-mhdx_tools-directory>  # Path to the mhdx_tools directory
"run_name": {library}_{pH}  # A descriptive name for your run
"names_and_seqs": {path-to-csv, input file 1}  # Path to the CSV input file
"timepoints": [timepoints in seconds]  # List of timepoints for your experiment
"runtime": <each's run length in minutes>  # Total runtime for each run
"time_bins": <a divisor of runtime>  # Time bin size dividing the runtime
```

There is also a block for file paths. Even if you transfered .gz files, here you should not use this suffix.
e.g. `20220212_Lib01_pH6_0s_01.mzML` not `20220212_Lib01_pH6_0s_01.mzML.gz`

A number of additional parameters is available in config/config.yaml. Default parameters were found empirally to perform best in our experiments. 

# **Running main mHDX-MS pipeline**

Our pipeline consists of three sequential steps. It was designed to take advantage of cluster resources. Edit the string after `--cluster` to match the cluster scheduler syntax from your cluster. Or omit if using local machine.

### Step 0: Activate python env `mhdxms`
```bash
conda activate mhdx
```

### **Step 1: Run preprocessing**
```bash
snakemake -s workflow/snakefiles/1_preprocessing.snakefile -j 1000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 04:00:00' --max-jobs-per-second 5
```

### **Step 2: Run factorization on undeuterated samples, check idotp and filter based on estimated qvalues**
```bash
snakemake -s workflow/snakefiles/2_idotp_check.snakefile -j 1000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 04:00:00' --max-jobs-per-second 5
```

### **Step 3: Run factorization and path optimizer on all timepoints**
```bash
snakemake -s workflow/snakefiles/3_main.snakefile -j 1000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 04:00:00' --max-jobs-per-second 5
```

### **Some of the expected final results**
```bash
resources/7_idotp_check/checked_library_info.json: Summary of identification results per protein.
resources/10_time_series_ics/consolidated_results.json: Path optimizer results.
results/plots/ic_time_series/ajf_plots: PDFs visualizing isotopic clusters.
results/plots/ic_time_series/winner_plots: PDFs of the final selected isotopic clusters.
```

## Authors

Állan Ferrari, Sugyan Dixit, Robert Wes Ludwig, Gabriel Rocklin

## Troubleshooting & Support ##
Ensure paths in config.yaml are correct before running.<br />
Adjust computational resources (`--mem=4GB, -t 04:00:00`) to match your HPC setup.<br />
For troubleshooting, review Snakemake logs and intermediate outputs.<br />
**Contact**: For any questions, contact `ajrferrari@gmail.com` or open an issue in this repository.<br />


# Demo

To facilitate reproducibility and provide a fully worked example, we offer a curated demo dataset at:

**[2025_Ferrari_mHDX_demo](https://nuwildcat.sharepoint.com/:f:/s/FSM-RocklinLab/EoPlHsYMMGNIlSeO1mWAl5ABpvMWeh918G7M7t0KBhES9A?e=IaLgqB)**

The demo dataset is organized as follows:

### Description of folders:

- **pH6/** and **pH9/**:
  - Contain all necessary processed files for each condition.
  - Includes isotope files (`.isotopes`), gzipped mzML files (`.mzML.gz`), and prepopulated `config.yaml` files.

- **sequences/**:
  - Contains `2018_HX_Mix1.csv`, a file with protein names, sequences, and monoisotopic masses.
  - This file should be used for setting up both the pH6 and pH9 experiments.

- **raw_files/**:
  - Contains the original raw mass spectrometry data.
  - Provided for users who wish to reprocess the raw data, generate their own mzML files using ProteoWizard, or recreate isotope files using IMTBX and Grppr.

### Important Note:

Before running the pipeline on the demo data, please edit the provided `config.yaml` files to update the `"mhdx_tools_dir"` field with the absolute path where you installed `mhdx_tools`.

All other settings are preconfigured and should work directly following the step-by-step instructions detailed above.

## System configuration

The pipeline was executed on the Northwestern Quest cluster equipped with the following specifications:

- **CPU**: Intel(R) Xeon(R) Gold 6230R CPU @ 2.10GHz
- **Memory**: 188 GB RAM per node
- **Operating System**: Fedora 8.10

## Pipeline benchmark

- **Total CPU time**: 52 hours
- **Execution mode**: Single-core execution (benchmark simulated)

