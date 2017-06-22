import errno
import os
import sys
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool

import dask
from dask.delayed import delayed
from numpy import asarray, empty, float64
from pandas import DataFrame
from tqdm import tqdm

from ._ffi import ffi
from ._ffi.lib import (close_bgen, close_variant_genotype, free, get_ncombs,
                       get_nsamples, get_nvariants, open_bgen,
                       open_variant_genotype, read_samples,
                       read_variant_genotype, read_variants,
                       sample_ids_presence, string_duplicate)

dask.set_options(pool=ThreadPool(cpu_count()))

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

PY3 = sys.version_info >= (3, )


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

        ncombs = get_ncombs(vg)
        g = empty((nsamples, ncombs), dtype=float64)

        pg = ffi.cast("real *", g.ctypes.data)
        read_variant_genotype(self._indexing[0], vg, pg)

        close_variant_genotype(self._indexing[0], vg)

        return g


def _read_genotype(indexing, nsamples, nvariants, nalleless, verbose):

    genotype = []
    rgv = ReadGenotypeVariant(indexing)

    for i in tqdm(range(nvariants), desc='variants', disable=not verbose):
        x = delayed(rgv)(nsamples, nalleless[i], i)
        genotype += [x]

    return asarray(genotype)


def read_bgen(filepath, verbose=True):
    r"""Read a given BGEN file.

    Args
    ----
    filepath : str
        A BGEN file path.
    verbose : bool
        ``True`` to show progress; ``False`` otherwise.

    Returns
    -------
    dict
        variants : Variant position, chromossomes, RSIDs, etc.
        samples : Sample identifications.
        genotype : Array of genotype references.
    """

    if PY3:
        try:
            filepath = filepath.encode()
        except AttributeError:
            pass

    if (not os.path.exists(filepath)):
        raise FileNotFoundError(errno.ENOENT,
                                os.strerror(errno.ENOENT), filepath)

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

    genotype = _read_genotype(indexing, nsamples, nvariants, nalleless,
                              verbose)

    return dict(variants=variants, samples=samples, genotype=genotype)
