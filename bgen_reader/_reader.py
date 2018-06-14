from os.path import join, dirname, basename
from ._misc import (
    make_sure_bytes,
    create_string,
    check_file_exist,
    check_file_readable,
)

import xarray as xr
import dask.array as da
from dask.delayed import delayed
from numpy import float64, nan, full, inf, asarray, concatenate, stack
from pandas import DataFrame, concat
from tqdm import tqdm
from ._metadata import try_read_variants_metadata_file

try:
    from functools import lru_cache
except ImportError:
    from ._pylru import lrudecorator as lru_cache

from ._ffi import ffi
from ._ffi.lib import (
    bgen_close,
    bgen_close_variant_genotype,
    bgen_free_samples,
    bgen_free_variants_metadata,
    bgen_ncombs,
    bgen_phased,
    bgen_missing,
    bgen_ploidy,
    bgen_nsamples,
    bgen_nvariants,
    bgen_open,
    bgen_open_variant_genotype,
    bgen_read_samples,
    bgen_read_variant_genotype,
    bgen_read_variants_metadata,
    bgen_sample_ids_presence,
    bgen_load_variants_metadata,
)


def _read_variants_from_bgen_file(bfile, index, v):
    variants = bgen_read_variants_metadata(bfile, index, v)
    if variants == ffi.NULL:
        raise RuntimeError("Could not read variants metadata.")
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
        variants = try_read_variants_metadata_file(bfile, mfile, index, v)
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


@lru_cache(2)
def _genotype_block(indexing, nsamples, variant_idx, nvariants):

    max_ncombs = -inf
    ncombss = []
    variants = []
    phased = []
    ploidy = []
    missing = []

    for i in range(variant_idx, variant_idx + nvariants):
        vg = bgen_open_variant_genotype(indexing[0], i)

        ncombs = bgen_ncombs(vg)
        max_ncombs = max(max_ncombs, ncombs)
        ncombss.append(ncombs)

        phased.append([bgen_phased(vg)] * nsamples)

        ploidy.append([bgen_ploidy(vg, j) for j in range(nsamples)])

        missing.append([bgen_missing(vg, j) for j in range(nsamples)])

        g = full((nsamples, ncombs), nan, dtype=float64)

        pg = ffi.cast("double *", g.ctypes.data)
        bgen_read_variant_genotype(indexing[0], vg, pg)

        bgen_close_variant_genotype(indexing[0], vg)

        variants.append(g)

    G = full((nvariants, nsamples, max_ncombs), nan, dtype=float64)

    for i in range(0, nvariants):
        G[i, :, :ncombss[i]] = variants[i]

    phased = asarray(phased, int)
    ploidy = asarray(ploidy, int)
    missing = asarray(missing, int)

    variant_idxs = range(variant_idx, variant_idx + nvariants)

    data = stack([phased, ploidy, missing], axis=2)

    coords = {
        "variant": variant_idxs,
        "sample": range(nsamples),
        "data": ["phased", "ploidy", "missing"],
    }
    dims = ("variant", "sample", "data")
    X = xr.DataArray(data, coords=coords, dims=dims)

    return G, X


def _read_genotype(indexing, nsamples, nvariants, nalleless, size, verbose):

    genotype = []
    X = []

    c = int((1024 * 1024 * size / 8) // nsamples)
    step = min(c, nvariants)
    tqdm_kwds = dict(desc="Variant mapping", disable=not verbose)

    kws = {"pure": True, "traverse": False}
    Gcall = delayed(lambda *args: _genotype_block(*args)[0], **kws)
    Xcall = delayed(lambda *args: _genotype_block(*args)[1], **kws)
    for i in tqdm(range(0, nvariants, step), **tqdm_kwds):
        size = min(step, nvariants - i)
        tup = indexing, nsamples, i, size
        shape0 = size, nsamples, nan
        shape1 = size, nsamples, 3

        g = da.from_delayed(Gcall(*tup), shape0, float64)
        genotype.append(g)

        x = da.from_delayed(Xcall(*tup), shape1, float64)
        X.append(x)

    a = da.concatenate(genotype, allow_unknown_chunksizes=True)
    b = da.concatenate(X, allow_unknown_chunksizes=True)
    return a, b


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

    G, X = _read_genotype(index, nsamples, nvariants, nalls, size, verbose)

    return dict(variants=variants, samples=samples, genotype=G, X=X)
