"""
BGEN file format reader
=======================

`BGEN <http://www.well.ox.ac.uk/~gav/bgen_format/>`_ is a file format for storing large
genetic datasets. It supports both unphased genotypes and phased haplotype data with
variable ploidy and number of alleles. It was designed to provides a compact data
representation without sacrificing variant access performance.

Functions
---------
create_metadata_file  Create variants metadata file.
read_bgen             Read a given BGEN file.
test                  Verify this package's integrity.

Documentation can be found at <https://github.com/limix/bgen-reader-py>.
"""

from __future__ import absolute_import

from ._dosage import (
    allele_expectation,
    convert_to_dosage,
    compute_dosage,
    allele_frequency,
)
from ._example import example_files
from ._metadata import create_metadata_file
from ._reader import read_bgen
from ._testit import test

try:
    from ._ffi import ffi as _
except Exception as e:
    msg = "\nIt is likely caused by a broken installation of this package."
    msg += "\nPlease, make sure you have a C compiler and try to uninstall"
    msg += "\nand reinstall the package again."
    e.msg = e.msg + msg
    raise e


__version__ = "2.0.8"

__all__ = [
    "__version__",
    "test",
    "read_bgen",
    "create_metadata_file",
    "convert_to_dosage",
    "allele_expectation",
    "example_files",
    "compute_dosage",
    "allele_frequency",
]
