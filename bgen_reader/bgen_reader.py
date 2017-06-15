# pylint: disable=E0401
from ._ffi.lib import (reader_close, reader_nsamples, reader_nvariants,
                       reader_open)


def print_something(filepath):
    bgenfile = reader_open(filepath)
    print("Number of samples: %d" % reader_nsamples(bgenfile))
    print("Number of variants: %d" % reader_nvariants(bgenfile))
    reader_close(bgenfile)
