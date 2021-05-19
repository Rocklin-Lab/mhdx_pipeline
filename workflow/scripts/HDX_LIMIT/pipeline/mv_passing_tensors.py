def main(input_paths, output_paths):
	for fin, fout in zip(input_paths, output_paths):
		os.makedirs(os.path.dirname(fout), exist_ok=True)
		shutil.copy(fin, fout)

if __name__ == "__main__":

	if "snakemake" in globals():
		main(snakemake.input, snakemake.output)

	else:
		parser = argparse.ArgumentParser()
		parser.add_argument("library_info_path", help="path/to/checked_library_info.json")
		parser.add_argument("input_dir_path", help="path/to/input_dir, inputs are globbed from this path, I don't work right now", default="resources/")
		parser.add_argument("output_dir_path", help="path/to/output_dir, outputs are globbed from this path, I don't work right now")
		args = parser.parse_args()

		# Make explicit inputs and outputs from checked_library_info.json.
		