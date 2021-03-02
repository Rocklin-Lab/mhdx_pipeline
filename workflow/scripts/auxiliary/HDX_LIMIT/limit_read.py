def limit_read(path):
    return cpickle.loads(zlib.decompress(open(path, "rb").read()))