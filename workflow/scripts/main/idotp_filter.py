import argparse
import pandas as pd

def main(all_idotp_csv_inputs, outpath, cutoff=0.95):
	out = []
	for fn in all_idotp_csv_inputs:
		lib_idx = int(fn.split('/')[-1].split('_')[0])
		idpc = pd.read_csv(fn)
		if idpc["idotp"].values[0] >= cutoff:
			out.append(lib_idx)
	out=sorted(out)
	outdf = pd.DataFrame.from_dict({"passing_indices"}: out)
	outdf.to_csv(outpath)
	return outdf

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="makes a .csv of all library_info indices with idotp >= cutoff, default 0.95")
	
	parser.add_argument("all_idotp_csv_inputs", help="list of ")
	parser.add_argument("outpath", help="path/to/filter_passing_indices.csv")
	parser.add_argument("-c", "--cutoff", type=float, default=0.95, help="lower limit on dot-product between theoretical integrated m/z of POI and int. m/z of observed signal in question. Float in range [0,1], default 0.95 ")

	args = parser.parse_arguments()

	main(args.all_idotp_csv_inputs, args.outpath, args.cutoff)