import dask.dataframe as dd
import pytest
from dask.delayed import Delayed
from numpy.testing import assert_, assert_allclose
from pandas import Series

from bgen_reader import (
    allele_expectation,
    allele_frequency,
    compute_dosage,
    example_files,
    read_bgen,
)


def test_read_bgem_interface():
    with example_files("haplotypes.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        assert_(isinstance(bgen, dict))
        assert_(isinstance(bgen["variants"], dd.DataFrame))
        assert_(isinstance(bgen["samples"], Series))
        assert_(isinstance(bgen["genotype"], list))
        assert_(isinstance(bgen["genotype"][0], Delayed))


def test_allele_expectation_interface():

    with example_files("haplotypes.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        with pytest.raises(ValueError):
            allele_expectation(bgen, 1)

    with example_files("complex.23bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        e = allele_expectation(bgen, 3)
        assert_allclose(
            e, [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 2.0, 0.0]]
        )


def test_allele_frequency_interface():

    with example_files("complex.23bits.bgen") as filepath:
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
    with example_files("complex.23bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        e = allele_expectation(bgen, 3)
        assert_allclose(compute_dosage(e), [0, 0, 0, 0])
        assert_allclose(compute_dosage(e, 0), [1.0, 2.0, 1.0, 0.0])
