from pathlib import Path

from cbgen import bgen_file
from pandas import Series, read_csv


def get_samples(bgen_filepath, verbose: bool) -> Series:
    with bgen_file(bgen_filepath) as bgen:

        if bgen.contain_samples:
            samples = bgen.read_samples()

        else:
            if verbose:
                print(
                    "Sample IDs are not present in this file."
                    "I will generate them on my own:"
                    " sample_1, sample_2, and so on."
                )
            samples = generate_samples(bgen.nsamples)

    return samples


def read_samples_file(sample_filepath: Path, verbose: bool) -> Series:
    if verbose:
        print(f"Sample IDs are read from {sample_filepath}.")

    samples = read_csv(sample_filepath, sep=" ", skiprows=[1]).iloc[:, 0].astype("str")
    return Series(samples, dtype=str, name="id")


def generate_samples(nsamples: int) -> Series:
    return Series([f"sample_{i}" for i in range(nsamples)], dtype=str, name="id")
