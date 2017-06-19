# pylint: disable=E0401
from numpy import int64

from ._ffi import ffi
from ._ffi.lib import (close_bgen, free, get_nsamples, get_nvariants,
                       open_bgen, read_samples, string_duplicate)


def _to_string(v):
    v = string_duplicate(v)
    return ffi.string(v.str, v.len).decode()


#
# def _read_variants(bgenfile):
#     from pandas import DataFrame
#
#     nvariants = reader_nvariants(bgenfile)
#
#     ids = ffi.new("string *[%d]" % nvariants)
#     rsids = ffi.new("string *[%d]" % nvariants)
#     chroms = ffi.new("string *[%d]" % nvariants)
#
#     positions = ffi.new("inti[%d]" % nvariants)
#     nalleless = ffi.new("inti[%d]" % nvariants)
#
#     reader_read_variants(bgenfile, ids, rsids, chroms, positions, nalleless)
#
#     data = dict(id=[], rsid=[], chrom=[], pos=[], nalleles=[])
#     for i in reversed(range(nvariants)):
#         data['id'].append(_to_string(ids[i]))
#         data['rsid'].append(_to_string(rsids[i]))
#         data['chrom'].append(_to_string(chroms[i]))
#
#         data['pos'].append(positions[i])
#         data['nalleles'].append(nalleless[i])
#
#     return DataFrame(data=data)
#
#
def _read_samples(bgenfile):
    from pandas import DataFrame

    nsamples = get_nsamples(bgenfile)
    samples = read_samples(bgenfile)

    py_ids = []
    for i in range(nsamples):
        py_ids.append(_to_string(samples[i]))

    return DataFrame(data=dict(id=py_ids))


#
# def _read_genotype_chunk(bgenfile, variant_start, variant_end):
#
#     nsamples = reader_nsamples(bgenfile)
#     X = zeros((nsamples, variant_end - variant_start), int64)
#
#     return X

# def _read_genotype(bgenfile):
#     import dask.array as da
#     from dask.delayed import delayed
#
#     nsamples = reader_nsamples(bgenfile)
#     nvariants = reader_nvariants(bgenfile)
#
#
#     variant_start = 0
#     variant_chunk = 10
#     genotype = []
#     while (variant_start < nvariants):
#         variant_end = min(variant_start + variant_chunk, nvariants)
#
#         x = delayed(_read_genotype_chunk)(bgenfile, variant_start, variant_end)
#
#         shape = (nsamples, variant_end - variant_start)
#
#         genotype += [da.from_delayed(x, shape, int64)]
#         variant_start = variant_end
#
#     genotype = da.concatenate(genotype, axis=1)
#
#
#     return genotype


def read(filepath):

    bgenfile = open_bgen(filepath)

    samples = _read_samples(bgenfile)
    # variants = _read_variants(bgenfile)

    # genotype = _read_genotype(bgenfile)
    genotype = None

    close_bgen(bgenfile)

    # return (variants, samples, genotype)
    return (None, samples, None)
