# pylint: disable=E0401
import dask.array as da
from dask.delayed import delayed
from numpy import int64, float64, empty
from pandas import DataFrame

from scipy.special import binom

from ._ffi import ffi
from ._ffi.lib import (close_bgen, free, get_nsamples, get_nvariants,
                       open_bgen, read_samples, read_variant_genotypes,
                       read_variants, sample_ids_presence, string_duplicate)


def _to_string(v):
    v = string_duplicate(v)
    return ffi.string(v.str, v.len).decode()


def _read_variants(bgenfile):
    indexing = ffi.new("VariantIndexing *[1]")
    nvariants = get_nvariants(bgenfile)
    variants = read_variants(bgenfile, indexing)

    data = dict(id=[], rsid=[], chrom=[], pos=[], nalleles=[])
    for i in range(nvariants):
        data['id'].append(_to_string(variants[i].id))
        data['rsid'].append(_to_string(variants[i].rsid))
        data['chrom'].append(_to_string(variants[i].chrom))

        data['pos'].append(variants[i].position)
        data['nalleles'].append(variants[i].nalleles)

    return (DataFrame(data=data), indexing)


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


def _read_genotype_variant(indexing, nsamples, nalleles, variant_idx):

    g = read_variant_genotype(indexing[0], nsamples, variant_idx)

    # ncombs = int(binom(nalleles + vg[0].ploidy - 1, nalleles - 1))
    # shape = (nsamples, ncombs)
    # g = empty(shape, dtype=float64)

    # g = vg[0].probabilities

    return g


def _read_genotype(indexing, nsamples, nvariants, nalleless):

    genotype = []
    import pdb; pdb.set_trace()
    g = read_variant_genotype(indexing[0], nsamples, 0)
    for i in range(nvariants):

        x = delayed(_read_genotype_variant)(indexing, nsamples, nalleless[i], i)

        genotype += [x]

    return genotype


def read(filepath):

    bgenfile = open_bgen(filepath)

    if sample_ids_presence(bgenfile) == 0:
        print("Sample ids are not present")
        samples = _generate_samples(bgenfile)
    else:
        samples = _read_samples(bgenfile)

    variants, indexing = _read_variants(bgenfile)
    nalleless = variants['nalleles'].values

    nsamples = samples.shape[0]
    nvariants = variants.shape[0]
    close_bgen(bgenfile)

    # genotype = _read_genotype(indexing, nsamples, nvariants, nalleless)
    genotype = None

    return (variants, samples, genotype)
