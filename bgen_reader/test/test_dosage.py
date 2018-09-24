from __future__ import unicode_literals

from bgen_reader import read_bgen, allele_expectation, example_files

from numpy import isnan, all
from numpy.testing import assert_, assert_equal, assert_allclose


def test_dosage_example_32bits():
    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)

        r = allele_expectation(bgen["genotype"], 2, 2)

        assert_allclose(r[5, 7], [1.9556273911044997, 0.044372608895500334])
        assert_(all(isnan(r[0, 0])))
        assert_equal(r.shape, (199, 500, 2))

        assert_allclose(
            allele_expectation(bgen["genotype"][5, 7], 2, 2),
            [1.9556273911044997, 0.044372608895500334],
        )

        assert_allclose(
            allele_expectation(bgen["genotype"][5], 2, 2)[7],
            [1.9556273911044997, 0.044372608895500334],
        )
