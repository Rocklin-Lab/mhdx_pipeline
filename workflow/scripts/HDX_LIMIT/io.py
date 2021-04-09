import zlib
import _pickle as cpickle

def limit_write(obj, outpath):
	""" Writes a Python object as a zlib-compressed pickle.

	Args:
		obj (any Python Object):
		outpath (string): path/to/file.cpickle.zlib to write out

	Returns:
		None

	"""
	with open(outpath, "wb") as file:
		file.write(zlib.compress(cpickle.dumps(obj)))

def limit_read(path):
	""" Reads a Python object from a zlib-compressed pickle.

	Args:
		path (string): path/to/file.cpickle.zlib to read

	Returns:
		(obj): some python object

	"""
	return cpickle.loads(zlib.decompress(open(path, "rb").read()))