### version 0.21.357

# The hdx_limit pipeline with Snakemake

## Installing hdx_limit

Step 1: Install anaconda. 

Download the linux installer from this link: https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh <br />
Installation instructions: https://docs.anaconda.com/anaconda/install/linux/

Step 2: Configure snakemake environment. Skip this step if already configured.

`conda install -n base -c conda-forge mamba` <br /> 
`conda activate base` <br /> 
`mamba create -c conda-forge -c bioconda -n snakemake snakemake=6.7.0`<br /> 

Step 3: Activate snamake environment and clone repository

`source activate snakemake` <br />
`git clone https://github.com/Rocklin-Lab/hdx_limit-pipeline.git ` <br /> 
`cd hdx_limit-pipeline`<br />
`git submodule init` <br /> 
`git submodule update` <br /> 

## Configuration

The pipeline expects three input files:

1) a csv file contaning protein name, protein sequence, and monoisotopic mass. This file should be stored at resources/0_names_seqs_masses <br />
e.g. <br />
`name,sequence,mono_mass` <br />
`PDB2PJV,HMAVGIGALFLGFLGAAGSTVGAASGGGKKKKK,3085.722267` <br />
`PDB2N92,HMSWLSKTAKKLENSAKKRISEGIAIAIQGGPR,3604.9987800000013` <br />
`PDB2LZP,HMDTEIIGGLTIPPVVALVVMSRFGFFAHLLPR,3632.972752` <br />

2) mzML files from HX-MS data. These files should be stored at resources/0_mzml. Alternatively, the user can provide already gziped versions of the mzML files. In this case, mzML.gz files should be stored at 2_mzml_gz.

3) isotope file for each HX-MS 0 timepoint. This file is generated from IMTBX. More information can be found at https://dmtavt.com/IMTBX.


4) edit config/config.yaml to contain the specific paths to those files along with the set of timepoints for your experiment. 

A number of additional parameters are available in config/config.yaml. Default parameters were found to best perform in our experiments. 

## Running hdx_limit with Snakemake

Our pipeline consists of three sequential steps. Edit the string after `--cluster` to match the cluster scheduler syntax from your cluster.

`snakemake -s workflow/snakefiles/1_preprocessing.snakefile --use-conda -j 4000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 03:00:00' --max-jobs-per-second 5`
 
`snakemake -s workflow/snakefiles/2_idotp_check.snakefile --use-conda -j 4000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 03:00:00' --max-jobs-per-second 5`

`snakemake -s workflow/snakefiles/3_main.snakefile --use-conda -j 4000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 03:00:00' --max-jobs-per-second 5`

## Authors

## Credits + Aknowledgments

## Sources

## License
