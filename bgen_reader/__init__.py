"""
*******************
bgen_reader package
*******************

Bgen file format reader.

"""

from __future__ import absolute_import as _

try:
    from ._ffi import ffi as _
except Exception as e:
    msg = "\nIt is likely caused by a broken installation of this package."
    msg += "\nPlease, make sure you have a C compiler and try to uninstall"
    msg += "\nand reinstall the package again."
    e.msg = e.msg + msg
    raise e

from ._reader import read_bgen
from ._metadata import create_metadata_file
from ._testit import test
from ._dosage import convert_to_dosage

__version__ = "2.0.3"

__all__ = [
    '__version__', 'test', 'read_bgen', 'create_metadata_file',
    'convert_to_dosage'
]
