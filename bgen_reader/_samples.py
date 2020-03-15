from pathlib import Path

from pandas import Series, read_csv

from ._bgen_file import bgen_file


def get_samples(bgen_filepath, verbose: bool) -> Series:
    with bgen_file(bgen_filepath) as bgen:

        if bgen.contain_samples:
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


def read_samples_file(sample_filepath: Path, verbose: bool):
    if verbose:
        print(f"Sample IDs are read from {sample_filepath}.")

    samples = read_csv(sample_filepath, sep=" ", skiprows=[1]).iloc[:, 0].astype("str")
    return Series(samples, dtype=str, name="id")


def _generate_sample_ids(nsamples: int):
    return Series([f"sample_{i}" for i in range(nsamples)], dtype=str, name="id")
