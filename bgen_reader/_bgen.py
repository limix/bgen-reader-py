from contextlib import contextmanager

from ._ffi import ffi, lib
from ._string import make_sure_bytes


@contextmanager
def bgen_file(filepath):
    bgen = lib.bgen_file_open(make_sure_bytes(filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {filepath}.")
    try:
        yield bgen
    finally:
        lib.bgen_file_close(bgen)


@contextmanager
def bgen_metafile(filepath):
    metafile = lib.bgen_metafile_open(make_sure_bytes(filepath))
    if metafile == ffi.NULL:
        raise RuntimeError(f"Could not open {filepath}.")
    try:
        yield metafile
    finally:
        lib.bgen_metafile_close(metafile)
