# pylint: disable=E0401
import dask
import dask.array as da
from dask.delayed import delayed
from numpy import int64, float64, empty, asarray
from pandas import DataFrame
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool

dask.set_options(pool=ThreadPool(cpu_count()))

from scipy.special import binom

from ._ffi import ffi
from ._ffi.lib import (close_bgen, free, get_nsamples, get_nvariants,
                       open_bgen, read_samples, open_variant_genotype,
                       variant_genotype_ncombs, close_variant_genotype,
                       read_variants, sample_ids_presence, string_duplicate,
                       read_variant_genotype)

from tqdm import tqdm

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


class ReadGenotypeVariant(object):
    def __init__(self, indexing):
        self._indexing = indexing

    def __call__(self, nsamples, nalleles, variant_idx):

        vg = open_variant_genotype(self._indexing[0], variant_idx)

        ncombs = variant_genotype_ncombs(vg)
        g = empty((nsamples, ncombs), dtype=float64)

        pg = ffi.cast("real *", g.ctypes.data)
        read_variant_genotype(self._indexing[0], vg, pg)

        close_variant_genotype(self._indexing[0], vg)

        return g


def _read_genotype(indexing, nsamples, nvariants, nalleless):

    genotype = []
    rgv = ReadGenotypeVariant(indexing)

    for i in tqdm(range(nvariants), desc='variants'):
        x = delayed(rgv)(nsamples, nalleless[i], i)
        genotype += [x]

    return asarray(genotype)


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

    genotype = _read_genotype(indexing, nsamples, nvariants, nalleless)

    return (variants, samples, genotype)
