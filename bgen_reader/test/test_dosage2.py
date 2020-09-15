import pytest
from numpy import all, isnan, logical_not, zeros
from numpy.testing import assert_allclose, assert_equal

from bgen_reader import example_filepath, open_bgen


def test_dosage_example_32bits():
    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen:
        e = bgen.allele_expectation([5, 0])
        assert_allclose(e[7, 0, :], [1.9556273911044997, 0.044372608895500334])
        assert all(isnan(e[0, 1, :]))
        assert_equal(e.shape, (500, 2, 2))


def test_freq():
    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen:
        variant_index = bgen.rsids == "RSID_6"
        e = bgen.allele_expectation(variant_index)
        f = bgen.allele_frequency(e)
        assert_allclose(f[0, 0], 229.23103218810434)
        assert_allclose(f[0, 1], 270.7689678118956)


def test_dosage1():
    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen:
        variant_index = 3
        e = bgen.allele_expectation(variant_index)
        # Compute the dosage when considering the allele
        # in position 1 as the reference/alternative one.
        alt_allele_index = 1
        dosage = e[..., alt_allele_index]
        # Print the dosage of the first five samples only.
        # print(dosage[:5])
        assert_allclose(dosage[:2, 0], [1.9618530841455453, 0.009826655967586362])


def test_error():
    filepath = example_filepath("complex.bgen")
    with open_bgen(filepath, allow_complex=True, verbose=False) as bgen:
        with pytest.raises(ValueError):
            bgen.allele_expectation()  # some phased

        with pytest.raises(ValueError):
            # different #'s of alleles
            bgen.allele_expectation(logical_not(bgen.phased))
        with pytest.raises(ValueError):
            # nonconstant ploidy
            bgen.allele_expectation(logical_not(bgen.phased) * (bgen.nalleles == 2))
        e = bgen.allele_expectation(
            logical_not(bgen.phased) * (bgen.nalleles == 2),
            assume_constant_ploidy=False,
        )
        f = bgen.allele_frequency(e)
        assert_allclose(e[-1, -1, :], [1.0, 3.0])
        assert_allclose(f[-1, :], [5.0, 3.0])


def test_zero_width():
    filepath = example_filepath("complex.bgen")
    with open_bgen(filepath, allow_complex=True, verbose=False) as bgen:
        for assume_constant_ploidy in [False, True]:
            e = bgen.allele_expectation(
                [],
                assume_constant_ploidy=assume_constant_ploidy,
            )
            f = bgen.allele_frequency(e)
            assert e.shape == (bgen.nsamples, 0, bgen.nalleles[0])
            assert f.shape == (0, bgen.nalleles[0])

            good_variants = logical_not(bgen.phased) * (bgen.nalleles == 2)
            e = bgen.allele_expectation(
                ([], good_variants),
                assume_constant_ploidy=assume_constant_ploidy,
            )
            f = bgen.allele_frequency(e)
            assert e.shape == (0, sum(good_variants), bgen.nalleles[0])
            assert_equal(
                f, zeros((sum(good_variants), bgen.nalleles[0]))
            )  # We define the freq of something with no samples as 0

            e = bgen.allele_expectation(
                ([], []),
                assume_constant_ploidy=assume_constant_ploidy,
            )
            f = bgen.allele_frequency(e)
            assert e.shape == (0, 0, bgen.nalleles[0])
            assert f.shape == (0, bgen.nalleles[0])


def test_dosage2():
    import numpy as np
    import pandas as pd

    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen:
        variant_index = [3]
        assert bgen.ids[variant_index] == "SNPID_5"
        assert bgen.rsids[variant_index] == "RSID_5"
        probs, missing, ploidy = bgen.read(
            variant_index, return_missings=True, return_ploidies=True
        )
        assert not np.any(missing)
        assert np.all(ploidy == 2)
        df1 = pd.DataFrame(
            {
                "sample": bgen.samples,
                "0": probs[:, 0, 0],
                "1": probs[:, 0, 1],
                "2": probs[:, 0, 2],
            }
        )
        # print(df1)
        assert_allclose(df1.iloc[-1, -1], 0.015471935508649781)
        alleles_per_variant = [
            allele_ids.split(",") for allele_ids in bgen.allele_ids[variant_index]
        ]
        e = bgen.allele_expectation(variant_index)
        f = bgen.allele_frequency(e)
        df2 = pd.DataFrame(
            {
                "sample": bgen.samples,
                alleles_per_variant[0][0]: e[:, 0, 0],
                alleles_per_variant[0][1]: e[:, 0, 1],
            }
        )
        # print(df2)  # doctest: +NORMALIZE_WHITESPACE
        assert_allclose(df2.iloc[-1, -1], 1.0152583189809832)
        alt_index = f[0, :].argmin()
        alt = alleles_per_variant[0][alt_index]
        dosage = e[:, 0, alt_index]
        df4 = pd.DataFrame({"sample": bgen.samples, f"alt={alt}": dosage})
        assert_allclose(df4.iloc[-1, -1], 1.0152583189809832)


if __name__ == "__main__":
    pytest.main([__file__])
