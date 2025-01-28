### version 0.25.028

# The mHDX-MS pipeline with Snakemake

## Installing hdx_limit

Step 1: Install Miniconda (if you don't have yet): 

```
https://docs.anaconda.com/miniconda/install/#quick-command-line-install
```

Step 2: Create a new environment to support the pipeline:

```conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n mhdxms snakemake==7.26.0 python==3.9
```

Step 3: Activate mhdxms environment, clone pipeline and executable repositories, and install executable into current environment

```conda activate mhdxms
git clone https://github.com/Rocklin-Lab/hdx_limit.git 
python -m pip install ./hdx_limit
```


## Setting up a new experiment

For each mHDX-MS experiment, clone a new pipeline repository

`git clone https://github.com/Rocklin-Lab/hdx_limit-pipeline.git`

The pipeline expects a set of input files:

1) a csv file contaning protein name, protein sequence, and monoisotopic mass. This file should be stored at resources/0_names_seqs_masses <br />
e.g. names_and_seqs.csv <br />
```name,sequence,mono_mass
PDB2PJV,HMAVGIGALFLGFLGAAGSTVGAASGGGKKKKK,3085.722267
PDB2N92,HMSWLSKTAKKLENSAKKRISEGIAIAIQGGPR,3604.9987800000013
PDB2LZP,HMDTEIIGGLTIPPVVALVVMSRFGFFAHLLPR,3632.972752
```

2) mzML files from HX-MS experiment. These files should be stored at resources/0_mzml. Preferably, the user can provide already gziped versions of the mzML files. In this case, mzML.gz files should be stored at resources/2_mzml_gz. Current pipeline module to read .gz files was fully tested with `ms-convert` from `ProteomeWizard 3.0.21193-ccb3e0136`

Some options I suggest:
```
- Options
Binary encoding precision: 32-bit
Use zlib compression: True
Package in gzip: True
All other Options not selected

- Filters
Zero Samples: Remove --> Add
```

4) isotope file for each undeuterated sample. This file is generated with `IMTBX` and `Grppr`. More information can be found at `https://dmtavt.com/IMTBX`.

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
"hdx_limit_dir": <path-to-hdx_limit-directory>  # Path to the hdx_limit directory
"run_name": {library}_{pH}  # A descriptive name for your run
"names_and_seqs": {path-to-csv, input file 1}  # Path to the CSV input file
"timepoints": [timepoints in seconds]  # List of timepoints for your experiment
"runtime": <each's run length in minutes>  # Total runtime for each run
"time_bins": <a divisor of runtime>  # Time bin size dividing the runtime
```

There is also a block for file paths. Even if you transfered .gz files, here you should not use this suffix.
e.g. `20220212_Lib01_pH6_0s_01.mzML` not `20220212_Lib01_pH6_0s_01.mzML.gz`

A number of additional parameters is available in config/config.yaml. Default parameters were found to best perform in our experiments. 

## Running main mHDX-MS pipeline

Our pipeline consists of three sequential steps. It was designed to take advantage of cluster resources. Edit the string after `--cluster` to match the cluster scheduler syntax from your cluster. Or omit if using local machine.

```
# Step 0: Activate python env `mhdxms`
conda activate mhdxms

# Step 1: Run preprocessing
snakemake -s workflow/snakefiles/1_preprocessing.snakefile -j 1000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 04:00:00' --max-jobs-per-second 5

# Step 2: Run factorization on undeuterated samples and check idotp
snakemake -s workflow/snakefiles/2_idotp_check.snakefile -j 1000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 04:00:00' --max-jobs-per-second 5

# Step 3: Run factorization and path optimizer on all timepoints
snakemake -s workflow/snakefiles/3_main.snakefile -j 1000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 04:00:00' --max-jobs-per-second 5

```

## Authors

Állan Ferrari, Sugyan Dixit, Robert Wes Ludwig, Gabriel Rocklin

## Credits + Aknowledgments

## Sources

## License
