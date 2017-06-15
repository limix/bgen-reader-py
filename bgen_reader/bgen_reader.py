# pylint: disable=E0401
from ._ffi.lib import (read_variants, reader_close, reader_nsamples,
                       reader_nvariants, reader_open)


def print_something(filepath):
    bgenfile = reader_open(filepath)

    nsamples = reader_nsamples(bgenfile)
    nvariants = reader_nvariants(bgenfile)

    print("Number of samples: %d" % nsamples)
    print("Number of variants: %d" % nvariants)

    # ids = ffi.new("byte **", nvariants)
    # read_variantid_blocks

    reader_close(bgenfile)
