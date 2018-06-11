"""
*******************
bgen_reader package
*******************

Bgen file format reader.

"""

from __future__ import absolute_import as _

from ._reader import read_bgen
from ._metadata import create_metadata_file
from ._testit import test

__version__ = "2.0.0"

__all__ = ['__version__', 'test', 'read_bgen',
           'create_metadata_file']
