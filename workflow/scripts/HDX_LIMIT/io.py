import zlib
import _pickle as cpickle


def limit_write(obj, out_path):
    """ Writes a Python object as a zlib-compressed pickle.

	Args:
		obj (any Python Object):
		outpath (string): path/to/file.cpickle.zlib to write out

	Returns:
		None

	"""
    with open(out_path, "wb") as file:
        file.write(zlib.compress(cpickle.dumps(obj)))


def limit_read(path):
    """ Reads a Python object from a zlib-compressed pickle.

	Args:
		path (string): path/to/file.cpickle.zlib to read

	Returns:
		(obj): some python object

	"""
    return cpickle.loads(zlib.decompress(open(path, "rb").read()))


def optimize_paths_inputs(name, library_info, config):
    # Pass inputs as fxn of rt-group name wildcard. Creates analyze_tensor() input filenames in fixed pattern, input tensor names include library_info.index and rt-group avg elution time.
    name_inputs = []
    idxs = library_info.index[library_info["name"] == name].tolist()
    for key in config["timepoints"]:
        if len(config[key]) > 1:
            for file in config[key]:
                for idx in idxs:
                    name_inputs.append(
                        "resources/subtensor_ics/"
                        + str(idx)
                        + "_"
                        + file
                        + ".gz.cpickle.zlib"
                    )  # TODO: This may break when using .raw as first input, investigate
        else:
            file = config[key][0]
            for idx in idxs:
                name_inputs.append(
                    "resources/subtensor_ics/"
                    + str(idx)
                    + "_"
                    + file
                    + ".gz.cpickle.zlib"
                )

    return name_inputs


def idotp_check_inputs(i):
    # Writes inputs for idotp_check rule based on library_info index
    idx_inputs = []
    if len(config[0]) > 1:
        for file in config[0]:
            idx_inputs.append(
                "resources/tensors/" + str(i) + "_" + file + ".gz.cpickle.zlib"
            )
    else:
        file = config[0][0]
        idx_inputs.append(
            "resources/tensors/" + str(i) + "_" + file + ".gz.cpickle.zlib"
        )
    return idx_inputs
