import pytest
from numpy.testing import assert_allclose, assert_equal

from bgen_reader import example_files, read_bgen
from bgen_reader._example._files import can_run_with


@pytest.mark.slowtest
@pytest.mark.skipif(
    not can_run_with(["large.bgen"]), reason="cannot run with large.bgen file"
)
def test_large_file(capsys):
    with example_files("large.bgen") as filepath:
        data = read_bgen(filepath, verbose=False)
        variants = data["variants"]
        samples = data["samples"]
        genotype = data["genotype"]
        df = variants[
            (variants["chrom"] == "3")
            & (variants["pos"] > 250000)
            & (variants["pos"] < 290000)
        ].compute()
        assert_equal(len(df), 14)
        assert_equal(df.iloc[0]["id"], "")
        assert_equal(df.iloc[0]["rsid"], "Human_STR_911968")
        assert_equal(df.iloc[0]["nalleles"], 2)
        assert_equal(df.iloc[0]["allele_ids"], "<CNV>,CACAC")
        assert_equal(df.iloc[0]["vaddr"], 73115140)

        assert_equal(samples[382], "SAMEA2586986")
        assert_equal(df.index[0], 283970)
        g = genotype[df.index[0]].compute()
        assert_equal(g["probs"].shape, (415, 3))

        assert_allclose(
            g["probs"][410, :], [0.0, 0.4268863965819791, 0.5731136034180209]
        )
        assert_equal(g["phased"], False)
        assert_equal(g["missing"][3], False)
        assert_equal(g["ploidy"][3], 2)
