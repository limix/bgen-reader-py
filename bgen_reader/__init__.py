"""
BGEN file format reader
=======================

`BGEN <http://www.well.ox.ac.uk/~gav/bgen_format/>`_ is a file format for storing large
genetic datasets. It supports both unphased genotypes and phased haplotype data with
variable ploidy and number of alleles. It was designed to provides a compact data
representation without sacrificing variant access performance.

Functions
---------
allele_expectation  Compute the expectation of each allele.
allele_frequency    Compute allele frequency from its expectation.
compute_dosage      Compute dosage from allele expectation.
create_metafile     Create metafile.
example_filepath    Get file path to a file example.
read_bgen           Read a given BGEN file (original Dask-inspired API).
open_bgen           Read a given BGEN file (new NumPy-inspired API).
test                Verify this package's integrity.

Documentation can be found at <https://github.com/limix/bgen-reader-py>.
"""
from ._bgen2 import open_bgen
from ._dosage import allele_expectation, allele_frequency, compute_dosage
from ._example import example_filepath
from ._metafile import create_metafile
from ._reader import read_bgen
from ._testit import test

__version__ = "4.0.8"

__all__ = [
    "__version__",
    "allele_expectation",
    "allele_frequency",
    "compute_dosage",
    "create_metafile",
    "example_filepath",
    "read_bgen",
    "open_bgen",
    "test",
]
