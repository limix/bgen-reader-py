import os
import sys
import warnings
from os.path import basename, dirname, join
from threading import RLock

import dask.array as da
import dask.dataframe as dd
import xarray as xr
from cachetools import LRUCache, cached
from dask.array import concatenate, from_delayed
from dask.delayed import delayed
from numpy import asarray, float64, full, inf, nan, stack
from pandas import DataFrame, Series, read_csv
from tqdm import tqdm

from ._ffi import ffi
from ._ffi.lib import (
    bgen_close,
    bgen_close_metafile,
    bgen_close_variant_genotype,
    bgen_contain_samples,
    bgen_free_samples,
    bgen_free_variants_metadata,
    bgen_load_variants_metadata,
    bgen_max_nalleles,
    bgen_metafile_nparts,
    bgen_missing,
    bgen_ncombs,
    bgen_nsamples,
    bgen_nvariants,
    bgen_open,
    bgen_open_metafile,
    bgen_open_variant_genotype,
    bgen_phased,
    bgen_ploidy,
    bgen_read_partition,
    bgen_read_samples,
    bgen_read_variant_genotype,
    bgen_read_variants_metadata,
)
from ._file import (
    _get_temp_filepath,
    assert_file_exist,
    assert_file_readable,
    permission_write_file,
)
from ._metadata import create_metafile
from ._misc import create_string, make_sure_bytes

PY3 = sys.version_info >= (3,)

if not PY3:
    FileNotFoundError = IOError


def read_bgen(filepath, metafile_filepath=None, samples_filepath=None, verbose=True):
    r"""Read a given BGEN file.

    Parameters
    ----------
    filepath : str
        A BGEN file path.
    size : float, optional
        Chunk size in megabytes. Defaults to ``50``.
    metafile_filepath : str, optional
        TODO: fix it
        If ``True``, it will try to read the variants metadata from the
        metadata file ``filepath + ".metadata"``. If this is not possible,
        the variants metadata will be read from the BGEN file itself. If
        ``filepath + ".metadata"`` does not exist, it will try to create one
        with the same name to speed up reads. If ``False``, variants metadata
        will be read only from the BGEN file. If a file path is given instead,
        it assumes that the specified metadata file is valid and readable and
        therefore it will read variants metadata from that file only. Defaults
        to ``True``.
    samples_filepath : str, optional
        A sample file in `GEN format <https://goo.gl/bCzo7m>`_.
        If samples_filepath is provided, sample IDs are read from this file. Otherwise,
        it reads from the BGEN file itself if present. Defaults to ``None``.
    verbose : bool, optional
        ``True`` to show progress; ``False`` otherwise.

    Returns
    -------
    variants : :class:`pandas.DataFrame`
        Variant position, chromossomes, RSIDs, etc.
    samples : :class:`pandas.DataFrame`
        Sample identifications.
    genotype : :class:`dask.array.Array`
        Array of genotype references.
    X : :class:`dask.array.Array`
        Allele probabilities.

    Note
    ----
    Metadata files can speed up subsequent reads tremendously. But often the user does
    not have write permission for the default metadata file location
    ``filepath + ".metadata"``. We thus provide the
    :function:`bgen_reader.create_metafile` function for creating one at the
    given path.
    """

    assert_file_exist(filepath)
    assert_file_readable(filepath)

    metafile_filepath = _get_valid_metafile_filepath(filepath, metafile_filepath)
    if not os.path.exists(metafile_filepath):
        create_metafile(filepath, metafile_filepath, verbose)

    samples = _get_samples(filepath, samples_filepath, verbose)
    variants = _map_metadata(filepath, metafile_filepath, samples)

    return dict(variants=variants, samples=samples)

    # nsamples = samples.shape[0]
    # bgen_close(bgen)

    # variants = _read_genotype2(index, nsamples, variants, verbose)
    # variants.name = "variants"

    # return dict(variants=variants, samples=samples)


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


def _read_variants2(bfile, filepath, metadata_file, verbose):
    if verbose:
        v = 1
    else:
        v = 0

    mfile = metadata_file
    mf = ffi.new("struct bgen_mf **")
    nvariants = bgen_nvariants(bfile)

    assert mfile
    if mfile is False:
        variants = _read_variants_from_bgen_file(bfile, mf, v)
    elif mfile is True:
        mfile = join(dirname(filepath), basename(filepath) + b".metadata")
        variants = try_read_variants_metadata_file2(bfile, mfile, mf, v)
    else:
        variants = bgen_load_variants_metadata(bfile, mfile, mf, v)

    if variants == ffi.NULL:
        msg = "Could not read the metadata file {}.".format(mfile)
        raise RuntimeError(msg)

    data = _create_variants_dataframe(variants, nvariants)
    bgen_free_variants_metadata(bfile, variants)

    return DataFrame(data=data), mf


def _get_npartitions(bgen_filepath, metafile_filepath):
    bgen = bgen_open(make_sure_bytes(bgen_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {bgen_filepath}.")

    metafile = bgen_open_metafile(make_sure_bytes(metafile_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {metafile_filepath}.")

    nparts = bgen_metafile_nparts(metafile)

    if bgen_close_metafile(metafile) != 0:
        raise RuntimeError(f"Error while closing metafile: {metafile_filepath}.")

    bgen_close(bgen)

    return nparts

def _get_nvariants(bgen_filepath):
    bgen = bgen_open(make_sure_bytes(bgen_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {bgen_filepath}.")

    nvariants = bgen_nvariants(bgen)
    bgen_close(bgen)

    return nvariants


def _map_metadata(bgen_filepath, metafile_filepath, samples):
    nparts = _get_npartitions(bgen_filepath, metafile_filepath)
    nvariants = _get_nvariants(bgen_filepath)
    dfs = []
    index_base = 0
    part_size = nvariants // nparts
    divisions = []
    for i in range(nparts):
        divisions.append(index_base)
        d = delayed(_read_partition)(bgen_filepath, metafile_filepath, i, index_base)
        dfs.append(d)
        index_base += part_size
    divisions.append(nvariants - 1)
    meta = [
        ("id", str),
        ("rsid", str),
        ("chrom", str),
        ("pos", int),
        ("nalleles", int),
        ("allele_ids", str),
    ]
    df = dd.from_delayed(dfs, meta=dd.utils.make_meta(meta), divisions=divisions)
    return df


def _bgen_str_to_str(s):
    if s.str == ffi.NULL:
        return ""
    return ffi.string(s.str, s.len).decode()


def _read_allele_ids(metadata):
    n = metadata.nalleles
    alleles = [_bgen_str_to_str(metadata.allele_ids[i]) for i in range(n)]
    return ",".join(alleles)


def _read_partition(bgen_filepath, metafile_filepath, part, index_base):
    bgen = bgen_open(make_sure_bytes(bgen_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {bgen_filepath}.")

    metafile = bgen_open_metafile(make_sure_bytes(metafile_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {metafile_filepath}.")

    nvariants_ptr = ffi.new("int *")
    metadata = bgen_read_partition(metafile, part, nvariants_ptr)
    nvariants = nvariants_ptr[0]
    variants = []
    for i in range(nvariants):
        id_ = _bgen_str_to_str(metadata[i].id)
        rsid = _bgen_str_to_str(metadata[i].rsid)
        chrom = _bgen_str_to_str(metadata[i].chrom)
        pos = metadata[i].position
        nalleles = metadata[i].nalleles
        allele_ids = _read_allele_ids(metadata[i])
        variants.append([id_, rsid, chrom, pos, nalleles, allele_ids])

    index = range(index_base, index_base + nvariants)
    variants = DataFrame(
        variants,
        index=index,
        columns=["id", "rsid", "chrom", "pos", "nalleles", "allele_ids"],
        dtype=str,
    )
    variants["pos"] = variants["pos"].astype(int)
    variants["nalleles"] = variants["nalleles"].astype(int)

    if bgen_close_metafile(metafile) != 0:
        raise RuntimeError(f"Error while closing metafile: {metafile_filepath}.")

    bgen_close(bgen)

    return variants


def _read_samples(bgen_file, verbose):
    if verbose:
        verbose = 1
    else:
        verbose = 0

    nsamples = bgen_nsamples(bgen_file)
    samples = bgen_read_samples(bgen_file, verbose)

    ids = [create_string(samples[i]) for i in range(nsamples)]

    bgen_free_samples(bgen_file, samples)
    return Series(ids, dtype=str, name="id")


def _read_samples_from_file(sample_file, verbose):
    if verbose:
        print("Sample IDs are read from {}.".format(sample_file))

    samples = read_csv(sample_file, sep=" ", skiprows=[1]).iloc[:, 0].astype("str")
    return Series(samples, dtype=str, name="id")


def _generate_samples(bgen_file):
    nsamples = bgen_nsamples(bgen_file)
    return Series([f"sample_{i}" for i in range(nsamples)], dtype=str, name="id")


cache = LRUCache(maxsize=2)
lock = RLock()


@cached(cache, lock=lock)
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
        G[i, :, : ncombss[i]] = variants[i]

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


# @cached(cache, lock=lock)
def _genotype_block2(indexing, nsamples, variant_idx):

    # open("/Users/horta/tmp/{}".format(variant_idx), "w").close()
    vg = bgen_open_variant_genotype(indexing[0], variant_idx)

    ncombs = bgen_ncombs(vg)
    phased = bgen_phased(vg)

    ploidy = [bgen_ploidy(vg, j) for j in range(nsamples)]
    missing = [bgen_missing(vg, j) for j in range(nsamples)]

    prob = full((nsamples, ncombs), nan, dtype=float64)

    bgen_read_variant_genotype(indexing[0], vg, ffi.cast("double *", prob.ctypes.data))
    bgen_close_variant_genotype(indexing[0], vg)

    prob = xr.DataArray(
        prob,
        dims=["sample", "probability"],
        coords={"ploidy": ("sample", ploidy), "missing": ("sample", missing)},
        attrs={"phased": phased},
    )

    return prob


def _read_genotype(indexing, nsamples, nvariants, nalleless, size, verbose):

    genotype = []
    X = []

    c = max(int((1024 * 1024 * size / 8) // nsamples), 1)
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


class VariantGenotype(object):
    def __init__(self, delayed):
        self._delayed = delayed

    def compute(self):
        return self._delayed.compute()

    def __repr__(self):
        return "<call compute()>"

    def __str__(self):
        return "<compute()>"


def _read_genotype2(indexing, nsamples, variants, verbose):
    tqdm_kwds = dict(desc="Variant mapping", disable=not verbose)
    genotype = []

    kws = {"pure": True, "traverse": False}
    shape = (len(variants),)
    call = delayed(lambda *args: _genotype_block2(*args), **kws)
    for i in tqdm(range(len(variants)), **tqdm_kwds):
        data = da.from_delayed(call(indexing, nsamples, i), shape, object)
        genotype.append(VariantGenotype(data))

    variants["genotype"] = Series(genotype)

    return variants


# def read_bgen_old(
#     filepath, size=50, verbose=True, metadata_file=True, sample_file=None
# ):
#     r"""Read a given BGEN file.

#     Parameters
#     ----------
#     filepath : str
#         A BGEN file path.
#     size : float, optional
#         Chunk size in megabytes. Defaults to ``50``.
#     verbose : bool, optional
#         ``True`` to show progress; ``False`` otherwise.
#     metadata_file : bool, str, optional
#         If ``True``, it will try to read the variants metadata from the
#         metadata file ``filepath + ".metadata"``. If this is not possible,
#         the variants metadata will be read from the BGEN file itself. If
#         ``filepath + ".metadata"`` does not exist, it will try to create one
#         with the same name to speed up reads. If ``False``, variants metadata
#         will be read only from the BGEN file. If a file path is given instead,
#         it assumes that the specified metadata file is valid and readable and
#         therefore it will read variants metadata from that file only. Defaults
#         to ``True``.
#     sample_file : str, optional
#         A sample file in `GEN format <https://goo.gl/bCzo7m>`_.
#         If sample_file is provided, sample IDs are read from this file. Otherwise, it
#         reads from the BGEN file itself if present. Defaults to ``None``.

#     Returns
#     -------
#     variants : :class:`pandas.DataFrame`
#         Variant position, chromossomes, RSIDs, etc.
#     samples : :class:`pandas.DataFrame`
#         Sample identifications.
#     genotype : :class:`dask.array.Array`
#         Array of genotype references.
#     X : :class:`dask.array.Array`
#         Allele probabilities.

#     Note
#     ----
#     Metadata files can speed up subsequent reads tremendously. But often the user does
#     not have write permission for the default metadata file location
#     ``filepath + ".metadata"``. We thus provide the
#     :function:`bgen_reader.create_metafile` function for creating one at the
#     given path.
#     """

#     filepath = make_sure_bytes(filepath)
#     if sample_file is not None:
#         sample_file = make_sure_str(sample_file)

#     assert_file_exist(filepath)
#     assert_file_readable(filepath)

#     if metadata_file not in [True, False]:
#         metadata_file = make_sure_bytes(metadata_file)
#         try:
#             assert_file_exist(metadata_file)
#         except FileNotFoundError as e:
#             msg = (
#                 "\n\nMetadata file `{}` does not exist.\nIf you want to create a "
#                 "metadata file in a custom location, please use "
#                 "`bgen_reader.create_metafile`.\n"
#             )
#             print(msg.format(metadata_file))
#             raise e
#         assert_file_readable(metadata_file)

#     bfile = bgen_open(filepath)
#     if bfile == ffi.NULL:
#         raise RuntimeError("Could not read {}.".format(filepath))

#     if sample_file is not None:
#         assert_file_exist(sample_file)
#         samples = _read_samples_from_file(sample_file, verbose)
#     elif bgen_sample_ids_presence(bfile) == 0:
#         if verbose:
#             print("Sample IDs are not present in this file.")
#             msg = "I will generate them on my own:"
#             msg += " sample_1, sample_2, and so on."
#             print(msg)
#         samples = _generate_samples(bfile)
#     else:
#         samples = _read_samples(bfile, verbose)

#     variants, index = _read_variants(bfile, filepath, metadata_file, verbose)
#     nalls = variants["nalleles"].values

#     nsamples = samples.shape[0]
#     nvariants = variants.shape[0]
#     bgen_close(bfile)

#     G, X = _read_genotype(index, nsamples, nvariants, nalls, size, verbose)

#     return dict(variants=variants, samples=samples, genotype=G, X=X)


_metafile_not_found = """\
Metafile `{filepath}` does not exist. If you wish to create a metafile in a custom
location, please use `bgen_reader.create_metafile`.
"""

_metafile_nowrite_dir = """\
You don't have permission to write to `{filepath}`.
This might prevent speeding-up the reading process in future runs.
"""


def _get_valid_metafile_filepath(bgen_filepath, metafile_filepath):
    if metafile_filepath is None:
        metafile = _infer_metafile_filenames(bgen_filepath)
        if os.path.exists(metafile["filepath"]):
            try:
                assert_file_readable(metafile["filepath"])
                return metafile["filepath"]
            except RuntimeError as e:
                warnings.warn(str(e), UserWarning)
                return _get_temp_filepath(metafile["dir"], metafile["filename"])
        else:
            if permission_write_file(metafile["filepath"]):
                return metafile["filepath"]
            else:
                fp = metafile["filepath"]
                warnings.warn(_metafile_nowrite_dir.format(filepath=fp), UserWarning)
                return _get_temp_filepath(metafile["dir"], metafile["filename"])
    else:
        try:
            assert_file_exist(metafile_filepath)
        except FileNotFoundError as e:
            fp = metafile_filepath
            warnings.warn(_metafile_not_found.format(filepath=fp), UserWarning)
            raise e
    return metafile_filepath


def _infer_metafile_filenames(bgen_filepath):
    metafile_filepath = os.path.abspath(bgen_filepath + ".metadata")
    metafile_filename = os.path.basename(metafile_filepath)
    metafile_dir = os.path.dirname(metafile_filepath)
    return {
        "filepath": metafile_filepath,
        "dir": metafile_dir,
        "filename": metafile_filename,
    }


def _get_samples(bgen_filepath, samples_filepath, verbose):
    bgen = bgen_open(make_sure_bytes(bgen_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {bgen_filepath}.")

    if samples_filepath is not None:
        assert_file_exist(samples_filepath)
        assert_file_readable(samples_filepath)
        samples = _read_samples_from_file(samples_filepath, verbose)
    elif bgen_contain_samples(bgen) == 0:
        if verbose:
            print("Sample IDs are not present in this file.")
            msg = "I will generate them on my own:"
            msg += " sample_1, sample_2, and so on."
            print(msg)
        samples = _generate_samples(bgen)
    else:
        samples = _read_samples(bgen, verbose)

    bgen_close(bgen)
    return samples
