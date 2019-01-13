"""
BGEN file format reader
=======================

`BGEN <http://www.well.ox.ac.uk/~gav/bgen_format/>`_ is a file format for storing large
genetic datasets. It supports both unphased genotypes and phased haplotype data with
variable ploidy and number of alleles. It was designed to provides a compact data
representation without sacrificing variant access performance.

Functions
---------
create_metafile  Create variants metadata file.
read_bgen        Read a given BGEN file.
test             Verify this package's integrity.

Documentation can be found at <https://github.com/limix/bgen-reader-py>.
"""
from ._dosage import allele_expectation, allele_frequency, compute_dosage
from ._example import example_files
from ._metadata import create_metafile
from ._reader import read_bgen
from ._testit import test

_ffi_err = """
It is likely caused by a broken installation of this package.
Please, make sure you have a C compiler and try to uninstall
and reinstall the package again."""

try:
    from ._ffi import ffi as _
except Exception as e:
    e.msg = e.msg + _ffi_err
    raise e

__version__ = "3.0.0"

__all__ = [
    "__version__",
    "test",
    "read_bgen",
    "create_metafile",
    "allele_expectation",
    "example_files",
    "compute_dosage",
    "allele_frequency",
]
