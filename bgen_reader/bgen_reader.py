from os import access, W_OK
from os.path import join, dirname, basename, exists
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
from ._misc import (
    make_sure_bytes,
    create_string,
    check_file_exist,
    check_file_readable,
)

import dask
import dask.array as da
from dask.delayed import delayed
from numpy import empty, float64, zeros
from pandas import DataFrame
from tqdm import tqdm

from ._ffi import ffi
from ._ffi.lib import (
    bgen_close,
    bgen_close_variant_genotype,
    bgen_free_samples,
    bgen_free_variants_metadata,
    bgen_ncombs,
    bgen_nsamples,
    bgen_nvariants,
    bgen_open,
    bgen_open_variant_genotype,
    bgen_read_samples,
    bgen_read_variant_genotype,
    bgen_read_variants_metadata,
    bgen_sample_ids_presence,
    bgen_create_variants_metadata_file,
    bgen_load_variants_metadata,
    bgen_store_variants_metadata,
)

dask.set_options(pool=ThreadPool(cpu_count()))


def _read_variants_from_bgen_file(bfile, index, v):
    variants = bgen_read_variants_metadata(bfile, index, v)
    if variants == ffi.NULL:
        raise RuntimeError("Could not read variants metadata.")
    return variants


def _try_read_variants_metadata_file(bfile, mfilepath, index, v):
    if exists(mfilepath):
        variants = bgen_load_variants_metadata(bfile, mfilepath, index, v)
        if variants == ffi.NULL:
            if v == 1:
                msg = "Warning: could not read variants"
                msg += " metadata from {}.".format(mfilepath)
                print(msg)
            variants = bgen_read_variants_metadata(bfile, index, v)
    else:
        variants = bgen_read_variants_metadata(bfile, index, v)

    if variants == ffi.NULL:
        raise RuntimeError("Could not read variants metadata.")

    errmsg = "Warning: could not create"
    errmsg += " the metadata file {}.".format(mfilepath)

    if not exists(mfilepath):
        if access("/path/to/folder", W_OK):
            e = bgen_store_variants_metadata(
                bfile, variants, index[0], mfilepath
            )
            if e != 0 and v == 1:
                print(errmsg)
        elif v == 1:
            print(errmsg)
    return variants


def _create_variants_dataframe(variants, nvariants):
    data = dict(id=[], rsid=[], chrom=[], pos=[], nalleles=[], allele_ids=[])
    for i in range(nvariants):
        data["id"].append(create_string(variants[i].id))
        data["rsid"].append(create_string(variants[i].rsid))
        data["chrom"].append(create_string(variants[i].chrom))

        data["pos"].append(variants[i].position)
        nalleles = variants[i].nalleles
        data["nalleles"].append(nalleles)
        alleles = []
        for j in range(nalleles):
            alleles.append(create_string(variants[i].allele_ids[j]))
        data["allele_ids"].append(",".join(alleles))
    return data


def _read_variants(bfile, filepath, metadata_file, verbose):
    if verbose:
        v = 1
    else:
        v = 0

    mfile = metadata_file
    index = ffi.new("struct bgen_vi **")
    nvariants = bgen_nvariants(bfile)

    if mfile is False:
        variants = _read_variants_from_bgen_file(bfile, index, v)
    elif mfile is True:
        mfile = join(dirname(filepath), basename(filepath) + b".metadata")
        variants = _try_read_variants_metadata_file(bfile, mfile, index, v)
    else:
        variants = bgen_load_variants_metadata(bfile, mfile, index, v)

    if variants == ffi.NULL:
        msg = "Could not read the metadata file {}.".format(mfile)
        raise RuntimeError(msg)

    data = _create_variants_dataframe(variants, nvariants)
    bgen_free_variants_metadata(bfile, variants)

    return (DataFrame(data=data), index)


def _read_samples(bgen_file, verbose):
    if verbose:
        verbose = 1
    else:
        verbose = 0

    nsamples = bgen_nsamples(bgen_file)
    samples = bgen_read_samples(bgen_file, verbose)

    ids = [create_string(samples[i]) for i in range(nsamples)]

    bgen_free_samples(bgen_file, samples)
    return DataFrame(data=dict(id=ids))


def _generate_samples(bgen_file):
    nsamples = bgen_nsamples(bgen_file)
    return DataFrame(data=dict(id=["sample_%d" % i for i in range(nsamples)]))


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
            G[i, :, : ncombss[i]] = variants[i]

        return G


def _read_genotype(indexing, nsamples, nvariants, nalleless, size, verbose):

    genotype = []
    rgv = ReadGenotypeVariant(indexing)

    c = int((1024 * 1024 * size / 8) // nsamples)
    step = min(c, nvariants)
    tqdm_kwds = dict(desc="Variant mapping", disable=not verbose)

    for i in tqdm(range(0, nvariants, step), **tqdm_kwds):
        size = min(step, nvariants - i)
        tup = nsamples, nalleless[i : i + size], i, size
        delayed_kwds = dict(pure=True, traverse=False)
        g = delayed(rgv, **delayed_kwds)(*tup)
        # TODO: THIS IS A HACK
        ncombs = 3
        g = da.from_delayed(g, (size, nsamples, ncombs), float64)
        genotype.append(g)

    return da.concatenate(genotype)


def read_bgen(filepath, size=50, verbose=True, metadata_file=True):
    r"""Read a given BGEN file.

    Parameters
    ----------
    filepath : str
        A BGEN file path.
    size : float
        Chunk size in megabytes. Defaults to ``50``.
    verbose : bool
        ``True`` to show progress; ``False`` otherwise.
    metadata_file : bool, str
        If ``True``, it will try to read the variants metadata from the
        metadata file ``filepath + ".metadata"``. If this is not possible,
        the variants metadata will be read from the BGEN file itself. If
        ``filepath + ".metadata"`` does not exist, it will try to create one
        with the same name to speed up reads. If ``False``, variants metadata
        will be read only from the BGEN file. If a file path is given instead,
        it assumes that the specified metadata file is valid and readable and
        therefore it will read variants metadata from that file only. Defaults
        to ``True``.

    Returns
    -------
    dict
        variants : Variant position, chromossomes, RSIDs, etc.
        samples : Sample identifications.
        genotype : Array of genotype references.
    """

    filepath = make_sure_bytes(filepath)

    check_file_exist(filepath)
    check_file_readable(filepath)

    if metadata_file not in [True, False]:
        metadata_file = make_sure_bytes(metadata_file)
        check_file_exist(metadata_file)
        check_file_readable(metadata_file)

    bfile = bgen_open(filepath)
    if bfile == ffi.NULL:
        raise RuntimeError("Could not read {}.".format(filepath))

    if bgen_sample_ids_presence(bfile) == 0:
        if verbose:
            print("Sample IDs are not present in this file.")
            msg = "I will generate them on my own:"
            msg += " sample_1, sample_2, and so on."
            print(msg)
        samples = _generate_samples(bfile)
    else:
        samples = _read_samples(bfile, verbose)

    variants, index = _read_variants(bfile, filepath, metadata_file, verbose)
    nalls = variants["nalleles"].values

    nsamples = samples.shape[0]
    nvariants = variants.shape[0]
    bgen_close(bfile)

    genotype = _read_genotype(index, nsamples, nvariants, nalls, size, verbose)

    return dict(variants=variants, samples=samples, genotype=genotype)


def create_metadata_file(bgen_filepath, metadata_filepath, verbose=True):
    r"""Create variants metadata file.

    Variants metadata file helps speed up subsequent reads of the associated
    BGEN file.

    Parameters
    ----------
    bgen_filepath : str
        BGEN file path.
    metadata_file : bool, str
        Metadata file path.
    verbose : bool
        ``True`` to show progress; ``False`` otherwise.
    """
    if verbose:
        verbose = 1
    else:
        verbose = 0

    bgen_filepath = make_sure_bytes(bgen_filepath)
    metadata_filepath = make_sure_bytes(metadata_filepath)

    check_file_exist(bgen_filepath)
    check_file_readable(bgen_filepath)

    if exists(metadata_filepath):
        raise ValueError(
            "The file {} already exists.".format(metadata_filepath)
        )

    e = bgen_create_variants_metadata_file(
        bgen_filepath, metadata_filepath, verbose
    )

    if e != 0:
        raise RuntimeError("Error while creating metadata file: {}".format(e))


def convert_to_dosage(G):
    r"""Convert probabilities to dosage.

    Let :math:`\mathbf G` be a three-dimensional array for which
    :math:`G_{i, j, l}` is the probability of the `j`-th sample having the
    `l`-th genotype (or haplotype) for the `i`-th locus.
    This function will return a bi-dimensional array ``X`` such that
    :math:`X_{i, j}` is the dosage of the `j`-th sample for the `i`-th locus.

    Parameters
    ----------
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
