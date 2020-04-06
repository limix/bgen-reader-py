from numpy import all, isnan
from numpy.testing import assert_allclose, assert_equal

from bgen_reader import allele_expectation, example_filepath, read_bgen


def test_dosage_example_32bits():
    filepath = example_filepath("example.32bits.bgen")
    bgen = read_bgen(filepath, verbose=False)

    e = allele_expectation(bgen, 5)
    assert_allclose(e[7], [1.9556273911044997, 0.044372608895500334])

    e = allele_expectation(bgen, 0)
    assert all(isnan(e[0]))

    e = allele_expectation(bgen, 0)
    assert_equal(e.shape, (500, 2))
