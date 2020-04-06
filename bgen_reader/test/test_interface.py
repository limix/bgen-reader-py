import dask.dataframe as dd
import pytest
from dask.delayed import Delayed
from numpy.testing import assert_allclose
from pandas import Series

from bgen_reader import (
    allele_expectation,
    allele_frequency,
    compute_dosage,
    example_filepath,
    read_bgen,
)


def test_read_bgem_interface():
    filepath = example_filepath("haplotypes.bgen")
    bgen = read_bgen(filepath, verbose=False)
    assert isinstance(bgen, dict)
    assert isinstance(bgen["variants"], dd.DataFrame)
    assert isinstance(bgen["samples"], Series)
    assert isinstance(bgen["genotype"], list)
    assert isinstance(bgen["genotype"][0], Delayed)


def test_allele_expectation_interface():
    bgen = read_bgen(example_filepath("haplotypes.bgen"), verbose=False)
    with pytest.raises(ValueError):
        allele_expectation(bgen, 1)

    bgen = read_bgen(example_filepath("complex.23bits.bgen"), verbose=False)
    e = allele_expectation(bgen, 3)
    assert_allclose(
        e, [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 2.0, 0.0]]
    )


def test_allele_frequency_interface():
    filepath = example_filepath("complex.23bits.bgen")
    with pytest.raises(ValueError):
        bgen = read_bgen(filepath, verbose=False)
        allele_expectation(bgen, 1)

    bgen = read_bgen(filepath, verbose=False)
    expec = allele_expectation(bgen, 3)
    freq = allele_frequency(expec)
    assert_allclose(freq, [1.33333333333, 1.0, 0.0])

    freq = allele_frequency(
        [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 2.0, 0.0]]
    )
    assert_allclose(freq, [1.33333333333, 1.0, 0.0])

    with pytest.raises(ValueError):
        allele_frequency([2, 3, 1])


def test_dosage_interface():
    bgen = read_bgen(example_filepath("complex.23bits.bgen"), verbose=False)
    e = allele_expectation(bgen, 3)
    assert_allclose(compute_dosage(e), [0, 0, 0, 0])
    assert_allclose(compute_dosage(e, 0), [1.0, 2.0, 1.0, 0.0])
