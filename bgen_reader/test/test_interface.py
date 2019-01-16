import dask.dataframe as dd
import pytest
from dask.delayed import Delayed
from numpy.testing import assert_, assert_allclose
from pandas import Series

from bgen_reader import allele_expectation, allele_frequency, example_files, read_bgen


def test_read_bgem_interface():
    with example_files("haplotypes.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        assert_(isinstance(bgen, dict))
        assert_(isinstance(bgen["variants"], dd.DataFrame))
        assert_(isinstance(bgen["samples"], Series))
        assert_(isinstance(bgen["genotype"], list))
        assert_(isinstance(bgen["genotype"][0], Delayed))


def test_allele_expectation():

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
