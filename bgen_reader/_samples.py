from pandas import Series, read_csv

from ._bgen_file import bgen_file
from ._file import assert_file_exist2, assert_file_readable2


def get_samples(bgen_filepath, samples_filepath, verbose: bool) -> Series:
    with bgen_file(bgen_filepath) as bgen:

        if samples_filepath is not None:
            assert_file_exist2(samples_filepath)
            assert_file_readable2(samples_filepath)
            samples = _read_samples_from_file(samples_filepath, verbose)

        elif bgen.contain_samples:
            samples = bgen.read_samples(verbose)

        else:
            if verbose:
                print(
                    "Sample IDs are not present in this file."
                    "I will generate them on my own:"
                    " sample_1, sample_2, and so on."
                )
            samples = _generate_sample_ids(bgen.nsamples)

    return samples


def _read_samples_from_file(sample_file, verbose: bool):
    if verbose:
        print(f"Sample IDs are read from {sample_file}.")

    samples = read_csv(sample_file, sep=" ", skiprows=[1]).iloc[:, 0].astype("str")
    return Series(samples, dtype=str, name="id")


def _generate_sample_ids(nsamples: int):
    return Series([f"sample_{i}" for i in range(nsamples)], dtype=str, name="id")
