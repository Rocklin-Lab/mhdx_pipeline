# version 0.21.348

# The hdx_limit pipeline with Snakemake:

## Installing hdx_limit:

conda install -n base -c conda-forge mamba 
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake=6.7.0

git clone https://github.com/Rocklin-Lab/hdx_limit-pipeline.git
git submodule init
git submodule update

## Configuration:

## Running hdx_limit with Snakemake:

snakemake -s workflow/snakefiles/1_preprocessing.snakefile --use-conda -j 4000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 03:00:00' --max-jobs-per-second 5
 
snakemake -s workflow/snakefiles/2_idotp_check.snakefile --use-conda -j 4000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 03:00:00' --max-jobs-per-second 5

snakemake -s workflow/snakefiles/3_main.snakefile --use-conda -j 4000 --keep-going --cluster 'sbatch -A p31346 -p short -N 1 -n 1 --mem=3GB -t 03:00:00' --max-jobs-per-second 5

## Authors:

## Credits + Aknowledgments:

## Sources:

## License:
