from pathlib import Path
from threading import RLock
from typing import List

from cachetools import LRUCache, cached
from cbgen import bgen_file, bgen_metafile
from cbgen.typing import Genotype
from tqdm import trange


def create_genotypes(bgen: bgen_file, metafile_filepath, verbose):
    nvariants = bgen.nvariants

    desc = "Mapping genotypes"
    return [
        _get_read_genotype(bgen, metafile_filepath)(i, dask_key_name=str(i))
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
            variants = p.variants
            end = min(variants.size, (sub_part + 1) * spart_size)
            vaddrs = tuple(p.variants.offset[start:end].tolist())
            g: List[Genotype] = read_genotype_partition(bgen_filepath, vaddrs)
            gm = g[m]
            return {
                "probs": gm.probability,
                "phased": gm.phased,
                "ploidy": gm.ploidy,
                "missing": gm.missing,
            }

    name = "read_genotype-" + tokenize(bytes(metafile_filepath))
    return delayed(read_genotype, name, True, None, False)


cache = LRUCache(maxsize=3)
lock = RLock()


@cached(cache, lock=lock)
def read_genotype_partition(bgen_filepath: Path, offsets):
    with bgen_file(bgen_filepath) as bgen:
        return [bgen.read_genotype(offset) for offset in offsets]


def _estimate_best_nsub_parts(nsamples, part_size):
    # Assume ideal block size, `bs`: 16MB
    # Assume 16 bytes per genotype per sample, `vs`
    # ideal nvariants to read: iv = bs / (vs * nsamples)
    # We then use iv to figure out in how many parts a partition will be subdivided
    # Let part_size be the number of variants in a partition
    # nsub_parts = min(int(part_size / iv), 1)
    bs = 16 * 1024 * 1024
    vs = 16
    iv = bs / (vs * nsamples)
    return max(int(part_size / iv), 1)


def _ceildiv(a, b):
    return -(-a // b)
