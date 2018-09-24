from __future__ import unicode_literals

import os

from bgen_reader import (
    read_bgen,
    create_metadata_file,
    convert_to_dosage,
    allele_expectation,
)

from numpy import isnan, all
from numpy.testing import assert_, assert_equal, assert_allclose

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def test_dosage_example_32bits():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    filepath = os.path.join(folder, b"example.32bits.bgen")
    bgen = read_bgen(filepath, verbose=False)
    variants = bgen["variants"]
    samples = bgen["samples"]

    r = allele_expectation(bgen["genotype"], 2, 2)
    print(r[5, 7, 0])
    print(r[5, 7, 1])
    assert_allclose(r[5, 7], [1.9570311913166734, 1.0014038002121737])
    assert_(all(isnan(r[0, 0])))
    assert_equal(r.shape, (199, 500, 2))

    assert_allclose(
        allele_expectation(bgen["genotype"][5, 7], 2, 2),
        [1.9570311913166734, 1.0014038002121737],
    )

    assert_allclose(
        allele_expectation(bgen["genotype"][5], 2, 2)[7],
        [1.9570311913166734, 1.0014038002121737],
    )
