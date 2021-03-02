import os
import sys
import argparse
from pymzml.utils.utils import index_gzip
from pymzml.run import Reader


def main(args):
    """
    Create an indexed gzip mzML file from a plain mzML.
    Usage: python3 gzip_mzml.py <path/to/mzml> <path/to/output>
    """
    with open(args.inpath) as fin:
        fin.seek(0, 2)
        max_offset_len = fin.tell()
        max_spec_no = Reader(args.inpath).get_spectrum_count() + 10

    index_gzip(
        args.inpath, args.outpath, max_idx=max_spec_no, idx_len=len(str(max_offset_len))
    )

if __name__ == "__main__":

    # set expected command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', help='path to .mzML input file, resources/mzml/*.mzML in pipeline context')
    parser.add_argument('outpath', help='path to .mzML.gz output file')
    # parse given arguments
    args = parser.parse_args()
      
    main(args)
    


