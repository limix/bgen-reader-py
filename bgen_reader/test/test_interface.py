import dask.dataframe as dd
import pytest
from dask.delayed import Delayed
from numpy.testing import assert_, assert_allclose
from pandas import Series

from bgen_reader._test_files import get_filepath
from bgen_reader import (
    allele_expectation,
    allele_frequency,
    compute_dosage,
    read_bgen,
)


def test_read_bgem_interface():
    filepath = get_filepath("haplotypes.bgen")
    bgen = read_bgen(filepath, verbose=False)
    assert_(isinstance(bgen, dict))
    assert_(isinstance(bgen["variants"], dd.DataFrame))
    assert_(isinstance(bgen["samples"], Series))
    assert_(isinstance(bgen["genotype"], list))
    assert_(isinstance(bgen["genotype"][0], Delayed))


def test_allele_expectation_interface():
    bgen = read_bgen(get_filepath("haplotypes.bgen"), verbose=False)
    with pytest.raises(ValueError):
        allele_expectation(bgen, 1)

    bgen = read_bgen(get_filepath("complex.23bits.bgen"), verbose=False)
    e = allele_expectation(bgen, 3)
    assert_allclose(
        e, [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 2.0, 0.0]]
    )


def test_allele_frequency_interface():
    filepath = get_filepath("complex.23bits.bgen")
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
    bgen = read_bgen(get_filepath("complex.23bits.bgen"), verbose=False)
    e = allele_expectation(bgen, 3)
    assert_allclose(compute_dosage(e), [0, 0, 0, 0])
    assert_allclose(compute_dosage(e, 0), [1.0, 2.0, 1.0, 0.0])
