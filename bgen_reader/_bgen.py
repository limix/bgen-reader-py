from contextlib import contextmanager

from ._ffi import ffi, lib
from ._string import make_sure_bytes


@contextmanager
def bgen_file(filepath):
    bgen = lib.bgen_open(make_sure_bytes(filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {filepath}.")
    try:
        yield bgen
    finally:
        lib.bgen_close(bgen)


@contextmanager
def bgen_metafile(filepath):
    mf = lib.bgen_open_metafile(make_sure_bytes(filepath))
    if mf == ffi.NULL:
        raise RuntimeError(f"Could not open {filepath}.")
    try:
        yield mf
    finally:
        lib.bgen_close_metafile(mf)
