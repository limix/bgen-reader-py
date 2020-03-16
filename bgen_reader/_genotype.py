from threading import RLock

from cachetools import LRUCache, cached
from numpy import asarray, float64, full, nan
from tqdm import trange

from ._bgen_file import bgen_file
from ._bgen_metafile import bgen_metafile, read_partition
from ._ffi import ffi, lib


def create_genotypes(bgen: bgen_file, metafile_filepath, verbose):
    nvariants = bgen.nvariants

    rg = _get_read_genotype(bgen, metafile_filepath)

    desc = "Mapping genotypes"
    return [
        rg(i, dask_key_name=str(i))
        for i in trange(nvariants, desc=desc, disable=not verbose)
    ]


def _get_read_genotype(bgen: bgen_file, metafile_filepath):
    from dask import delayed
    from dask.base import tokenize

    nsamples = bgen.nsamples
    bgen_filepath = bgen.filepath

    def read_genotype(i: int):

        with bgen_metafile(metafile_filepath) as mf:
            part_size = mf.partition_size
            part = i // part_size
            j = i % part_size
            p = mf.read_partition(part)
            nsub_parts = _estimate_best_nsub_parts(nsamples, part_size)
            spart_size = max(1, part_size // nsub_parts)
            sub_part = j // spart_size
            m = j % spart_size
            start = sub_part * spart_size
            end = min(len(p), (sub_part + 1) * spart_size)
            vaddrs = tuple(p.iloc[start:end]["vaddr"].tolist())
            g = read_genotype_partition(bgen_filepath, vaddrs)
            return g[m]

    name = "read_genotype-" + tokenize(bytes(metafile_filepath))
    return delayed(read_genotype, name, True, None, False)


cache = LRUCache(maxsize=3)
lock = RLock()


@cached(cache, lock=lock)
def read_genotype_partition(bgen_filepath, vaddrs):
    genotypes = []
    for vaddr in vaddrs:
        with bgen_file(bgen_filepath) as bgen:
            genotype = bgen.read_genotype(vaddr)
            genotypes.append(genotype)
    return genotypes


def _estimate_best_nsub_parts(nsamples, part_size):
    # Assume ideal block size, `bs`: 256KB
    # Assume 16 bytes per genotype per sample, `vs`
    # ideal nvariants to read: iv = bs / (vs * nsamples)
    # We then use iv to figure out in how many parts a partition will be subdivided
    # Let part_size be the number of variants in a partition
    # nsub_parts = min(int(part_size / iv), 1)
    bs = 256 * 1024
    vs = 16
    iv = bs / (vs * nsamples)
    return max(int(part_size / iv), 1)


def _ceildiv(a, b):
    return -(-a // b)
