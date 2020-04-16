from numpy import all, isnan
from numpy.testing import assert_allclose, assert_equal
import pytest

from bgen_reader import example_filepath, open_bgen


def test_dosage_example_32bits():
    filepath = example_filepath("example.32bits.bgen")
    bgen = open_bgen(filepath, verbose=False)

    e = bgen.allele_expectation([5,0])
    assert_allclose(e[7,0,:], [1.9556273911044997, 0.044372608895500334])
    assert all(isnan(e[0,1,:]))
    assert_equal(e.shape, (500, 2, 2))

if __name__ == "__main__":  #!!!cmk99 remove?
    pytest.main([__file__])
