import argparse
import pandas as pd

def main(all_idotp_csv_inputs, outpath=None, return_flag=False, cutoff=0.95):
	"""
	Reads rt-group idotp

    Parameters:
    argument1 (int): Description of arg1

    Returns:
    int:Returning value

    """

	out_dict = {}

	filter_passing_indices = []
	for fn in all_idotp_csv_inputs:
		lib_idx = int(fn.split('/')[-1].split('_')[0])
		idpc = pd.read_csv(fn)
		if idpc["idotp"].values[0] >= cutoff:
			filter_passing_indices.append(lib_idx)

	# re-order indices
	filter_passing_indices=sorted(filter_passing_indices)
	# add passing indices to output dict
	out_dict["filter_passing_indices"] = filter_passing_indices
	# make df output option
	out_df = pd.DataFrame.from_dict(out_dict)
	
	if outpath is not None:
		out_df.to_csv(outpath)
	
	if return_flag:
		return out_dict

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="makes a .csv of all library_info indices with idotp >= cutoff, default 0.95")
	
	parser.add_argument("all_idotp_csv_inputs", help="list of ")
	parser.add_argument("-o", "--outpath", help="path/to/filter_passing_indices.csv")
	parser.add_argument("-c", "--cutoff", type=float, default=0.95, help="lower limit on dot-product between theoretical integrated m/z of POI and int. m/z of observed signal in question. Float in range [0,1], default 0.95 ")

	args = parser.parse_arguments()

	main(args.all_idotp_csv_inputs, outpath=args.outpath, cutoff=args.cutoff)