import os
import sys
import psutil
import pandas as pd

sys.path.append(os.getcwd() + "/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx

# open library_info
library_info = pd.read_csv(snakemake.input[0])

# Find timepoint of passed filename by config comparison
my_fn = snakemake.input[1]
for tp in snakemake.config["timepoints"]:
    for fn in snakemake.config[tp]:
        if fn in my_fn:
            my_tp = tp

process = psutil.Process(os.getpid())

# memory before init
print("Pre-Initialization: " + str(process.memory_info().rss / (1024 * 1024 * 1024)))

# init TG
tg = hx.TensorGenerator(my_fn, my_tp, library_info)

# profile memory after init
print("Post-Initialization: " + str(process.memory_info().rss / (1024 * 1024 * 1024)))

# factorize internal DT
tg.DataTensor.factorize(gauss_params=(3, 1))

# profile memory after factorization
print("Post-Factorization: " + str(process.memory_info().rss / (1024 * 1024 * 1024)))

# output ICs as flat list
all_ics = []
for factor in tg.DataTensor.factors:
    for ic in factor.isotope_clusters:
        all_ics.append(ic)

hx.limit_write(all_ics, snakemake.output[0])
