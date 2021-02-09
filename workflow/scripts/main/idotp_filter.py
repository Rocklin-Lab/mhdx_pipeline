import pandas as pd
filtered_library_info = pd.read_csv(snakemake.input.pop(0))
for fn in snakemake.input:
	lib_idx = int(fn.split('/')[-1].split('_')[0])
	idpc = pd.read_csv(fn)
	if idpc["idotp"].values[0] < snakemake.config["idotp_cutoff"]:
		filtered_library_info.drop(lib_idx)
filtered_library_info.to_csv(snakemake.output[0])