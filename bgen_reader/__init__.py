"""
*******************
bgen_reader package
*******************

Bgen file format reader.

"""

from __future__ import absolute_import as _absolute_import

from pkg_resources import DistributionNotFound as _DistributionNotFound
from pkg_resources import get_distribution as _get_distribution

from .bgen_reader import convert_to_dosage, read_bgen

try:
    __version__ = _get_distribution('bgen_reader').version
except _DistributionNotFound:
    __version__ = 'unknown'


def test():
    import os
    p = __import__('bgen_reader').__path__[0]
    src_path = os.path.abspath(p)
    old_path = os.getcwd()
    os.chdir(src_path)

    try:
        return_code = __import__('pytest').main(['-q', '--doctest-modules'])
    finally:
        os.chdir(old_path)

    if return_code == 0:
        print("Congratulations. All tests have passed!")

    return return_code


__all__ = ['test', 'read_bgen', 'convert_to_dosage']
