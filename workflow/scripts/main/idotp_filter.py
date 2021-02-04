print("0")
import copy
import ipdb
import pandas as pd
print("1")
filtered_library_info = pd.read_csv(snakemake.input.pop(0))
print("2")
ins = copy.copy(snakemake.input)
for fn in ins:
	lib_idx = int(fn.split('/')[-1].split('_')[0])
	try:
		if pd.read_csv(fn)["idotp"] < config["idotp_cutoff"]:
			filtered_library_info.drop(lib_idx)
	except:
		ipdb.set_trace()
filtered_library_info.to_csv(snakemake.output[0])