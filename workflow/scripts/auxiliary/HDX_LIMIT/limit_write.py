def limit_write(obj, outpath):
    with open(outpath, "wb") as file:
        file.write(zlib.compress(cpickle.dumps(obj)))