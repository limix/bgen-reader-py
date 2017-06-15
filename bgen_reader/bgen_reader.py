# pylint: disable=E0401
from ._ffi import ffi
from ._ffi.lib import (free, reader_close, reader_nsamples, reader_nvariants,
                       reader_open, reader_read_variants, reader_read_samples)


def _to_string(v):
    v = ffi.gc(v, free)
    return ffi.string(v[0].s, v[0].len).decode()


def _read_variants(bgenfile):
    from pandas import DataFrame

    nvariants = reader_nvariants(bgenfile)

    ids = ffi.new("string *[%d]" % nvariants)
    rsids = ffi.new("string *[%d]" % nvariants)
    chroms = ffi.new("string *[%d]" % nvariants)

    positions = ffi.new("inti[%d]" % nvariants)
    nalleless = ffi.new("inti[%d]" % nvariants)

    reader_read_variants(bgenfile, ids, rsids, chroms, positions, nalleless)

    data = dict(id=[], rsid=[], chrom=[], pos=[], nalleles=[])
    for i in reversed(range(nvariants)):
        data['id'].append(_to_string(ids[i]))
        data['rsid'].append(_to_string(rsids[i]))
        data['chrom'].append(_to_string(chroms[i]))

        data['pos'].append(positions[i])
        data['nalleles'].append(nalleless[i])

    return DataFrame(data=data)

def _read_samples(bgenfile):
    from pandas import DataFrame

    nsamples = reader_nsamples(bgenfile)

    ids = ffi.new("string *[%d]" % nsamples)

    reader_read_samples(bgenfile, ids)

    py_ids = []
    for i in range(nsamples):
        py_ids.append(_to_string(ids[i]))

    return DataFrame(data=dict(id=py_ids))


def read(filepath):

    bgenfile = reader_open(filepath)

    variants = _read_variants(bgenfile)
    samples = _read_samples(bgenfile)

    reader_close(bgenfile)

    return (variants, samples, None)
