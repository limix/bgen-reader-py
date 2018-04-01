import errno
import os
import stat
import sys
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool

import dask
import dask.array as da
from dask.delayed import delayed
from numpy import empty, float64, zeros
from pandas import DataFrame
from tqdm import tqdm

from ._ffi import ffi
from ._ffi.lib import (bgen_close, bgen_close_variant_genotype,
                       bgen_free_samples, bgen_free_variants_metadata,
                       bgen_ncombs, bgen_nsamples, bgen_nvariants, bgen_open,
                       bgen_open_variant_genotype, bgen_read_samples,
                       bgen_read_variant_genotype, bgen_read_variants_metadata,
                       bgen_sample_ids_presence)

dask.set_options(pool=ThreadPool(cpu_count()))

PY3 = sys.version_info >= (3, )

if not PY3:
    FileNotFoundError = IOError


def _create_string(v):
    s = ffi.new("char[]", v.len)
    ffi.memmove(s, v.str, v.len)
    return ffi.string(s, v.len).decode()


def _read_variants(bgen_file):
    verbose = 0
    indexing = ffi.new("struct bgen_vi **")
    nvariants = bgen_nvariants(bgen_file)
    variants = bgen_read_variants_metadata(bgen_file, indexing, verbose)

    data = dict(id=[], rsid=[], chrom=[], pos=[], nalleles=[], allele_ids=[])
    for i in range(nvariants):
        data['id'].append(_create_string(variants[i].id))
        data['rsid'].append(_create_string(variants[i].rsid))
        data['chrom'].append(_create_string(variants[i].chrom))

        data['pos'].append(variants[i].position)
        nalleles = variants[i].nalleles
        data['nalleles'].append(nalleles)
        alleles = []
        for j in range(nalleles):
            alleles.append(_create_string(variants[i].allele_ids[j]))
        data['allele_ids'].append(','.join(alleles))

    bgen_free_variants_metadata(bgen_file, variants)

    return (DataFrame(data=data), indexing)


def _read_samples(bgen_file):

    verbose = 0
    nsamples = bgen_nsamples(bgen_file)
    samples = bgen_read_samples(bgen_file, verbose)

    py_ids = []
    for i in range(nsamples):
        py_ids.append(_create_string(samples[i]))

    bgen_free_samples(bgen_file, samples)
    return DataFrame(data=dict(id=py_ids))


def _generate_samples(bgen_file):
    nsamples = bgen_nsamples(bgen_file)
    return DataFrame(data=dict(id=['sample_%d' % i for i in range(nsamples)]))


class ReadGenotypeVariant(object):
    def __init__(self, indexing):
        self._indexing = indexing

    def __call__(self, nsamples, nalleles, variant_idx, nvariants):

        ncombss = []
        variants = []

        for i in range(variant_idx, variant_idx + nvariants):
            vg = bgen_open_variant_genotype(self._indexing[0], i)

            ncombs = bgen_ncombs(vg)
            ncombss.append(ncombs)
            g = empty((nsamples, ncombs), dtype=float64)

            pg = ffi.cast("double *", g.ctypes.data)
            bgen_read_variant_genotype(self._indexing[0], vg, pg)

            bgen_close_variant_genotype(self._indexing[0], vg)

            variants.append(g)

        G = zeros((nvariants, nsamples, max(ncombss)), dtype=float64)

        for i in range(0, nvariants):
            G[i, :, :ncombss[i]] = variants[i]

        return G


def _read_genotype(indexing, nsamples, nvariants, nalleless, size, verbose):

    genotype = []
    rgv = ReadGenotypeVariant(indexing)

    c = int((1024 * 1024 * size / 8) // nsamples)
    step = min(c, nvariants)
    tqdm_kwds = dict(desc='Variant mapping', disable=not verbose)

    for i in tqdm(range(0, nvariants, step), **tqdm_kwds):
        size = min(step, nvariants - i)
        tup = nsamples, nalleless[i:i + size], i, size
        delayed_kwds = dict(pure=True, traverse=False)
        g = delayed(rgv, **delayed_kwds)(*tup)
        # TODO: THIS IS A HACK
        ncombs = 3
        g = da.from_delayed(g, (size, nsamples, ncombs), float64)
        genotype.append(g)

    return da.concatenate(genotype)


def read_bgen(filepath, size=50, verbose=True):
    r"""Read a given BGEN file.

    Args
    ----
    filepath : str
        A BGEN file path.
    size : float
        Chunk size in megabytes. Defaults to ``50``.
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
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                filepath)

    if not _group_readable(filepath):
        msg = "You don't have file"
        msg += " permission for reading {}.".format(filepath)
        raise RuntimeError(msg)

    bgen_file = bgen_open(filepath)
    if bgen_file == ffi.NULL:
        raise RuntimeError("Could not read {}.".format(filepath))

    if bgen_sample_ids_presence(bgen_file) == 0:
        if verbose:
            print("Sample IDs are not present in this file.")
            msg = "I will generate them on my own:"
            msg += " sample_1, sample_2, and so on."
            print(msg)
        samples = _generate_samples(bgen_file)
    else:
        samples = _read_samples(bgen_file)

    sys.stdout.write("Reading variants (it should take less than a minute)...")
    sys.stdout.flush()
    variants, indexing = _read_variants(bgen_file)
    sys.stdout.write(" done.\n")
    sys.stdout.flush()
    nalleless = variants['nalleles'].values

    nsamples = samples.shape[0]
    nvariants = variants.shape[0]
    bgen_close(bgen_file)

    genotype = _read_genotype(indexing, nsamples, nvariants, nalleless, size,
                              verbose)

    return dict(variants=variants, samples=samples, genotype=genotype)


def convert_to_dosage(G, verbose=True):
    r"""Convert probabilities to dosage.

    Let :math:`\mathbf G` be a three-dimensional array for which
    :math:`G_{i, j, l}` is the probability of the `j`-th sample having the
    `l`-th genotype (or haplotype) for the `i`-th locus.
    This function will return a bi-dimensional array ``X`` such that
    :math:`X_{i, j}` is the dosage of the `j`-th sample for the `i`-th locus.

    Args
    ----
    G : array_like
        A three-dimensional array.

    Returns
    -------
    dask_array
        Matrix representing dosages.
    """
    ncombs = G.shape[2]
    mult = da.arange(ncombs, chunks=ncombs, dtype=float64)
    return da.sum(mult * G, axis=2)


def _group_readable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IRGRP)
