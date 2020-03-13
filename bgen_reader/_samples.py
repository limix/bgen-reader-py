from pandas import Series, read_csv

from ._ffi import lib, ffi
from ._file import assert_file_exist, assert_file_readable
from ._bgen import bgen_file
from ._string import create_string


def get_samples(bgen_filepath, samples_filepath, verbose):
    with bgen_file(bgen_filepath) as bgen:
        if samples_filepath is not None:
            assert_file_exist(samples_filepath)
            assert_file_readable(samples_filepath)
            samples = _read_samples_from_file(samples_filepath, verbose)
        elif lib.bgen_file_contain_samples(bgen):
            samples = _read_samples(bgen, verbose)
        else:
            if verbose:
                print("Sample IDs are not present in this file.")
                msg = "I will generate them on my own:"
                msg += " sample_1, sample_2, and so on."
                print(msg)
            samples = _generate_samples(bgen)

    return samples


def _read_samples(bgen, verbose):
    if verbose:
        verbose = 1
    else:
        verbose = 0

    nsamples = lib.bgen_file_nsamples(bgen)
    samples = lib.bgen_file_read_samples(bgen, verbose)
    if samples == ffi.NULL:
        raise RuntimeError(f"Could not fetch samples from the bgen file.")

    ids = [create_string(lib.bgen_samples_get(samples, i)) for i in range(nsamples)]

    lib.bgen_samples_destroy(samples)
    return Series(ids, dtype=str, name="id")


def _read_samples_from_file(sample_file, verbose):
    if verbose:
        print("Sample IDs are read from {}.".format(sample_file))

    samples = read_csv(sample_file, sep=" ", skiprows=[1]).iloc[:, 0].astype("str")
    return Series(samples, dtype=str, name="id")


def _generate_samples(bgen):
    nsamples = lib.bgen_file_nsamples(bgen)
    return Series([f"sample_{i}" for i in range(nsamples)], dtype=str, name="id")
