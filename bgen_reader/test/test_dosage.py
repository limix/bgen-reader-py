from __future__ import unicode_literals

from numpy import all, isnan
from numpy.testing import assert_, assert_allclose, assert_equal

from bgen_reader import allele_expectation, example_files, read_bgen


def test_dosage_example_32bits():
    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)

        e = allele_expectation(bgen, 5)
        assert_allclose(e[7], [1.9556273911044997, 0.044372608895500334])

        e = allele_expectation(bgen, 0)
        assert_(all(isnan(e[0])))

        e = allele_expectation(bgen, 0)
        assert_equal(e.shape, (500, 2))
