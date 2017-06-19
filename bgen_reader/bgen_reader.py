# pylint: disable=E0401
from numpy import int64
from pandas import DataFrame

from ._ffi import ffi
from ._ffi.lib import (close_bgen, free, get_nsamples, get_nvariants,
                       open_bgen, read_samples, read_variants,
                       sample_ids_presence, string_duplicate)


def _to_string(v):
    v = string_duplicate(v)
    return ffi.string(v.str, v.len).decode()


def _read_variants(bgenfile):
    indexing = ffi.new("VariantIndexing *[1]")
    nvariants = get_nvariants(bgenfile)
    variants = read_variants(bgenfile, indexing)

    data = dict(id=[], rsid=[], chrom=[], pos=[], nalleles=[])
    for i in reversed(range(nvariants)):
        data['id'].append(_to_string(variants[i].id))
        data['rsid'].append(_to_string(variants[i].rsid))
        data['chrom'].append(_to_string(variants[i].chrom))

        data['pos'].append(variants[i].position)
        data['nalleles'].append(variants[i].nalleles)

    return DataFrame(data=data)


def _read_samples(bgenfile):

    nsamples = get_nsamples(bgenfile)
    samples = read_samples(bgenfile)

    py_ids = []
    for i in range(nsamples):
        py_ids.append(_to_string(samples[i]))

    return DataFrame(data=dict(id=py_ids))


def _generate_samples(bgenfile):
    nsamples = get_nsamples(bgenfile)
    return DataFrame(data=dict(id=['sample_%d' % i for i in range(nsamples)]))


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

    if sample_ids_presence(bgenfile) == 0:
        print("Sample ids are not present")
        samples = _generate_samples(bgenfile)
    else:
        samples = _read_samples(bgenfile)

    variants = _read_variants(bgenfile)

    # genotype = _read_genotype(bgenfile)
    genotype = None

    close_bgen(bgenfile)

    return (variants, samples, None)
