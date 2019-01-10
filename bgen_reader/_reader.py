import os
import warnings

import dask.dataframe as dd
from dask.delayed import delayed

from ._ffi import ffi
from ._ffi.lib import bgen_close, bgen_nvariants, bgen_open, bgen_read_variants_metadata
from ._file import (
    _get_temp_filepath,
    assert_file_exist,
    assert_file_readable,
    permission_write_file,
)
from ._metadata import create_metafile
from ._misc import bgen_str_to_str, create_string, make_sure_bytes
from ._partition import get_npartitions, read_partition
from ._samples import get_samples


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

    samples = get_samples(filepath, samples_filepath, verbose)
    variants = _map_metadata(filepath, metafile_filepath, samples)

    return dict(variants=variants, samples=samples)


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


def _get_nvariants(bgen_filepath):
    bgen = bgen_open(make_sure_bytes(bgen_filepath))
    if bgen == ffi.NULL:
        raise RuntimeError(f"Could not open {bgen_filepath}.")

    nvariants = bgen_nvariants(bgen)
    bgen_close(bgen)

    return nvariants


def _map_metadata(bgen_filepath, metafile_filepath, samples):
    nparts = get_npartitions(bgen_filepath, metafile_filepath)
    nvariants = _get_nvariants(bgen_filepath)
    dfs = []
    index_base = 0
    part_size = nvariants // nparts
    divisions = []
    for i in range(nparts):
        divisions.append(index_base)
        d = delayed(read_partition)(bgen_filepath, metafile_filepath, i, index_base)
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


def _map_genotype():
    pass
    # genotype[5] = {


# "probs": array([]),
# "phased": 0,
# "ploidy": array([2, 1, 2]),
# "missing": array([0, 0, 0]),
# }




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
