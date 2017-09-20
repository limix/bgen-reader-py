"""
*******************
bgen_reader package
*******************

Bgen file format reader.

"""

from __future__ import absolute_import as _

import _cffi_backend as _

from ._test import test
from .bgen_reader import convert_to_dosage, read_bgen

__name__ = "bgen-reader"
__version__ = "1.0.0"
__author__ = "Danilo Horta"
__author_email__ = "horta@ebi.ac.uk"

__all__ = [
    "__name__", "__version__", "__author__", "__author_email__", "test",
    'read_bgen', 'convert_to_dosage'
]
