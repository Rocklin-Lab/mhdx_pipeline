import os
import sys
from pymzml.utils.utils import index_gzip
from pymzml.run import Reader


def main(mzml_path, out_path):
    """
    Create an indexed gzip mzML file from a plain mzML.
    Usage: python3 gzip_mzml.py <path/to/mzml> <path/to/output>
    """
    with open(mzml_path) as fin:
        fin.seek(0, 2)
        max_offset_len = fin.tell()
        max_spec_no = Reader(mzml_path).get_spectrum_count() + 10

    index_gzip(
        mzml_path, out_path, max_idx = max_spec_no, idx_len = len(str(max_offset_len))
    )

try:
    main(snakemake.input[0], snakemake.output[0])
except:
        print(main.__doc__)
