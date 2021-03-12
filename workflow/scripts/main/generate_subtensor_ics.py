import os
import sys
import psutil
import argparse
import pandas as pd

sys.path.append(os.getcwd() + "/workflow/scripts/auxiliary/")
import LC_IM_MS_TensorAnalysis as hx

def main(library_info_path, tensor_input_path, isotope_clusters_output_path, timepoints, gauss_params=(3,1)):
	# open library_info
	library_info = pd.read_csv(library_info_path)

	# Find timepoint of passed filename by config comparison
	for tp in timepoints["timepoints"]:
	    for fn in timepoints[tp]:
	        if fn in tensor_input_path:
	            my_tp = tp

	process = psutil.Process(os.getpid())

	# memory before init
	print("Pre-Initialization: " + str(process.memory_info().rss / (1024 * 1024 * 1024)))

	# init TG
	tg = hx.TensorGenerator(tensor_input_path, my_tp, library_info)

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

	hx.limit_write(all_ics, isotope_clusters_output_path)

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="Accepts tensor as input, factorizes and saves IsotopeClusters from resulting Factors")

	parser.add_argument("library_info_path", help="path/to/library_info.csv")
	parser.add_argument("tensor_input_path", help="path/to/file.cpickle.zlib for tensor to factorize")
	parser.add_argument("isotope_clusters_output_path", help="path/to/output.cpickle.zlib, list of IsotopeClusters")
	parser.add_argument("timepoints", help="dictionary with 'timepoints' containing hdx times in seconds, and a key for each timepoint corresponding to a list of timepoint mzml filenames. Can pass opened snakemake.config object")
	parser.add_argument("-g", "--gauss_params", type=tuple, default=(3,1), help="determines intensity of gaussian smoothing in rt, and dt dimensions")

	main(library_info_path=args.library_info_path, tensor_input_path=args.tensor_input_path, isotope_clusters_output_path=args.isotope_clusters_output_path, timepoints=args.timepoints, gauss_params=args.gauss_params)