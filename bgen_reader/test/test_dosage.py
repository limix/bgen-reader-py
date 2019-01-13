from __future__ import unicode_literals

from bgen_reader import read_bgen, allele_expectation, example_files

from numpy import isnan, all
from numpy.testing import assert_, assert_equal, assert_allclose


def test_dosage_example_32bits():
    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)

        e = allele_expectation(bgen, 5)
        assert_allclose(e[7], [1.9556273911044997, 0.044372608895500334])

        e = allele_expectation(bgen, 0)
        assert_(all(isnan(e[0])))

        e = allele_expectation(bgen, 0)
        assert_equal(e.shape, (500, 2))

