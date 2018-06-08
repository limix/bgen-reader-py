"""
*******************
bgen_reader package
*******************

Bgen file format reader.

"""

from __future__ import absolute_import as _

from .bgen_reader import convert_to_dosage, create_metadata_file, read_bgen
from .testit import test

__version__ = "1.1.4"

__all__ = ['__version__', 'test', 'read_bgen',
           'convert_to_dosage', 'create_metadata_file']
