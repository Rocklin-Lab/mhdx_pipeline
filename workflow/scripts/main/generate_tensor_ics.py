import os
import sys
import yaml
import psutil
import argparse
import pandas as pd

sys.path.append(os.getcwd() + "/workflow/scripts/auxiliary/")
from HDX_LIMIT import TensorGenerator, limit_write

def main(library_info_path, tensor_input_path, timepoints_dict, isotope_clusters_output_path=None, return_flag=False, gauss_params=(3,1)):
	out = {}

	# open library_info
	library_info = pd.read_csv(library_info_path)

	# Find timepoint of passed filename by config comparison
	for tp in timepoints_dict["timepoints"]:
	    for fn in timepoints_dict[tp]:
	        if fn in tensor_input_path:
	            my_tp = tp

	process = psutil.Process(os.getpid())

	# memory before init
	print("Pre-Initialization: " + str(process.memory_info().rss / (1024 * 1024 * 1024)))

	# init TG
	tg = TensorGenerator(tensor_input_path, my_tp, library_info)

	# profile memory after init
	print("Post-Initialization: " + str(process.memory_info().rss / (1024 * 1024 * 1024)))

	# factorize internal DT
	tg.DataTensor.factorize(gauss_params=gauss_params)

	# profile memory after factorization
	print("Post-Factorization: " + str(process.memory_info().rss / (1024 * 1024 * 1024)))

	# output ICs as flat list
	all_ics = []
	for factor in tg.DataTensor.factors:
	    for ic in factor.isotope_clusters:
	        all_ics.append(ic)

	out['all_ics'] = all_ics

	if isotope_clusters_output_path is not None: 
		limit_write(all_ics, isotope_clusters_output_path)

	if return_flag:
		return out

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="Accepts tensor as input, factorizes and saves IsotopeClusters from resulting Factors")

	parser.add_argument("library_info_path", help="path/to/library_info.csv")
	parser.add_argument("tensor_input_path", help="path/to/file.cpickle.zlib for tensor to factorize")
	parser.add_argument("timepoints_yaml", help="path/to/file.yaml containing list of hdx timepoints in integer seconds which are also keys mapping to lists of each timepoint's .mzML file, can pass config/config.yaml - for Snakemake context")
	parser.add_argument("-o", "--isotope_clusters_output_path", help="path/to/output.cpickle.zlib, list of IsotopeClusters")
	parser.add_argument("-r", "--return_flag", type=bool, default=False, help="option to return output dictionary in python, for notebook context")
	parser.add_argument("-g", "--gauss_params", type=tuple, default=(3,1), help="determines intensity of gaussian smoothing in rt, and dt dimensions")

	args = parser.parse_args()

	#open .yaml into dict for main()
	open_timepoints = yaml.load(open(args.timepoints_yaml, 'rb').read())

	main(library_info_path=args.library_info_path, tensor_input_path=args.tensor_input_path, timepoints_dict=open_timepoints, isotope_clusters_output_path=args.isotope_clusters_output_path, return_flag=args.return_flag, gauss_params=args.gauss_params)