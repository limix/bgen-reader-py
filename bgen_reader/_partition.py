from pandas import DataFrame

from ._ffi import ffi
from ._ffi.lib import (
    bgen_close_metafile,
    bgen_metafile_nparts,
    bgen_open_metafile,
    bgen_read_partition,
)
from ._file import bgen_file
from ._misc import bgen_str_to_str, make_sure_bytes


def read_partition(bgen_filepath, metafile_filepath, part, index_base):
    with bgen_file(bgen_filepath):

        metafile = bgen_open_metafile(make_sure_bytes(metafile_filepath))
        if metafile == ffi.NULL:
            raise RuntimeError(f"Could not open {metafile_filepath}.")

        nvariants_ptr = ffi.new("int *")
        metadata = bgen_read_partition(metafile, part, nvariants_ptr)
        nvariants = nvariants_ptr[0]
        variants = []
        for i in range(nvariants):
            id_ = bgen_str_to_str(metadata[i].id)
            rsid = bgen_str_to_str(metadata[i].rsid)
            chrom = bgen_str_to_str(metadata[i].chrom)
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

    return variants


def get_npartitions(bgen_filepath, metafile_filepath):
    with bgen_file(bgen_filepath) as bgen:
        metafile = bgen_open_metafile(make_sure_bytes(metafile_filepath))
        if bgen == ffi.NULL:
            raise RuntimeError(f"Could not open {metafile_filepath}.")

        nparts = bgen_metafile_nparts(metafile)

        if bgen_close_metafile(metafile) != 0:
            raise RuntimeError(f"Error while closing metafile: {metafile_filepath}.")

    return nparts


def _read_allele_ids(metadata):
    n = metadata.nalleles
    alleles = [bgen_str_to_str(metadata.allele_ids[i]) for i in range(n)]
    return ",".join(alleles)
