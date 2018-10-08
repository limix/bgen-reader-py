from __future__ import unicode_literals

import os
import sys

import pytest
from bgen_reader import read_bgen, create_metadata_file, example_files
from numpy import isnan
from numpy.testing import assert_, assert_equal, assert_allclose

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def test_bgen_reader_variants_info():
    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"]
        samples = bgen["samples"]
        assert_("genotype" in bgen)

        assert_equal(variants.loc[0, "chrom"], "01")
        assert_equal(variants.loc[0, "id"], "SNPID_2")
        assert_equal(variants.loc[0, "nalleles"], 2)
        assert_equal(variants.loc[0, "allele_ids"], "A,G")
        assert_equal(variants.loc[0, "pos"], 2000)
        assert_equal(variants.loc[0, "rsid"], "RSID_2")

        assert_equal(variants.loc[7, "chrom"], "01")
        assert_equal(variants.loc[7, "id"], "SNPID_9")
        assert_equal(variants.loc[7, "nalleles"], 2)
        assert_equal(variants.loc[7, "allele_ids"], "A,G")
        assert_equal(variants.loc[7, "pos"], 9000)
        assert_equal(variants.loc[7, "rsid"], "RSID_9")

        n = variants.shape[0]
        assert_equal(variants.loc[n - 1, "chrom"], "01")
        assert_equal(variants.loc[n - 1, "id"], "SNPID_200")
        assert_equal(variants.loc[n - 1, "nalleles"], 2)
        assert_equal(variants.loc[n - 1, "allele_ids"], "A,G")
        assert_equal(variants.loc[n - 1, "pos"], 100001)
        assert_equal(variants.loc[n - 1, "rsid"], "RSID_200")

        assert_equal(samples.loc[0, "id"], "sample_001")
        assert_equal(samples.loc[7, "id"], "sample_008")

        n = samples.shape[0]
        assert_equal(samples.loc[n - 1, "id"], "sample_500")

        G = bgen["genotype"].compute()

        assert_(all(isnan(G[0, 0, :])))
        a = [0.027802362811705648, 0.00863673794284387, 0.9635608992454505]
        assert_allclose(G[0, 1, :], a)
        b = [
            0.97970582847010945215516,
            0.01947019668749305418287,
            0.00082397484239749366197,
        ]
        assert_allclose(G[1, 2, :], b)

    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"]
        assert_("samples" in bgen)
        assert_("genotype" in bgen)

        assert_equal(variants.loc[0, "chrom"], "01")
        G = bgen["genotype"].compute()
        assert_(all(isnan(G[0, 0, :])))
        a = [0.027802362811705648, 0.00863673794284387, 0.9635608992454505]
        assert_allclose(G[0, 1, :], a)

    if os.path.exists(filepath + ".metadata"):
        os.remove(filepath + ".metadata")


def test_bgen_reader_phased_genotype():
    _test_bgen_reader_phased_genotype(50)
    _test_bgen_reader_phased_genotype(0.0001)


def _test_bgen_reader_phased_genotype(size):
    with example_files("haplotypes.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False, size=size)
        variants = bgen["variants"]
        samples = bgen["samples"]
        assert_("genotype" in bgen)

        assert_equal(variants.loc[0, "chrom"], "1")
        assert_equal(variants.loc[0, "id"], "SNP1")
        assert_equal(variants.loc[0, "nalleles"], 2)
        assert_equal(variants.loc[0, "allele_ids"], "A,G")
        assert_equal(variants.loc[0, "pos"], 1)
        assert_equal(variants.loc[0, "rsid"], "RS1")

        assert_equal(variants.loc[2, "chrom"], "1")
        assert_equal(variants.loc[2, "id"], "SNP3")
        assert_equal(variants.loc[2, "nalleles"], 2)
        assert_equal(variants.loc[2, "allele_ids"], "A,G")
        assert_equal(variants.loc[2, "pos"], 3)
        assert_equal(variants.loc[2, "rsid"], "RS3")

        assert_equal(samples.loc[0, "id"], "sample_0")
        assert_equal(samples.loc[2, "id"], "sample_2")

        n = samples.shape[0]
        assert_equal(samples.loc[n - 1, "id"], "sample_3")

        G = bgen["genotype"].compute()
        a = [1.0, 0.0, 1.0, 0.0]
        assert_allclose(G[0, 0, :], a)
        k = len(variants)
        n = len(samples)
        a = [1.0, 0.0, 0.0, 1.0]
        assert_allclose(G[k - 1, n - 1, :], a)

    with example_files("haplotypes.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"]
        assert_("samples" in bgen)
        assert_("genotype" in bgen)

        assert_equal(variants.loc[0, "chrom"], "1")
        G = bgen["genotype"].compute()
        a = [1.0, 0.0, 1.0, 0.0]
        assert_allclose(G[0, 0, :], a)
        k = len(variants)
        n = len(samples)
        a = [1.0, 0.0, 0.0, 1.0]
        assert_allclose(G[k - 1, n - 1, :], a)

    if os.path.exists(filepath + ".metadata"):
        os.remove(filepath + ".metadata")


def test_bgen_reader_without_metadata():
    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False, metadata_file=False)
        variants = bgen["variants"]
        samples = bgen["samples"]
        assert_("genotype" in bgen)
        assert_equal(variants.loc[7, "allele_ids"], "A,G")
        n = samples.shape[0]
        assert_equal(samples.loc[n - 1, "id"], "sample_500")
        assert_(not os.path.exists(filepath + ".metadata"))


def test_bgen_reader_with_wrong_metadata_file():
    with example_files(["example.32bits.bgen", "wrong.metadata"]) as filepaths:
        with pytest.raises(RuntimeError):
            read_bgen(filepaths[0], verbose=False, metadata_file=filepaths[1])


def test_bgen_reader_with_nonexistent_metadata_file():
    with example_files("example.32bits.bgen") as filepath:
        folder = os.path.dirname(filepath)
        metadata_file = os.path.join(folder, "nonexistent.metadata")

        with pytest.raises(FileNotFoundError):
            read_bgen(filepath, verbose=False, metadata_file=metadata_file)


def test_bgen_reader_file_notfound():
    with pytest.raises(FileNotFoundError):
        read_bgen("/1/2/3/example.32bits.bgen", verbose=False)


def test_create_metadata_file():
    with example_files("example.32bits.bgen") as filepath:
        folder = os.path.dirname(filepath)
        metadata_file = os.path.join(folder, filepath + ".metadata")

        create_metadata_file(filepath, metadata_file, verbose=False)
        assert_(os.path.exists(metadata_file))
        os.remove(metadata_file)


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires python3.6 or higher")
def test_create_metadata_file_bytes():
    with example_files("example.32bits.bgen") as filepath:
        folder = os.path.dirname(filepath)
        metadata_file = os.path.join(folder, filepath + ".metadata")

        create_metadata_file(filepath, metadata_file.encode(), verbose=False)
        assert_(os.path.exists(metadata_file))
        os.remove(metadata_file)


def test_bgen_reader_complex():
    with example_files("complex.23bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"]
        samples = bgen["samples"]
        assert_("genotype" in bgen)

        assert_equal(variants.loc[0, "chrom"], "01")
        assert_equal(variants.loc[0, "id"], "")
        assert_equal(variants.loc[0, "nalleles"], 2)
        assert_equal(variants.loc[0, "allele_ids"], "A,G")
        assert_equal(variants.loc[0, "pos"], 1)
        assert_equal(variants.loc[0, "rsid"], "V1")

        assert_equal(variants.loc[7, "chrom"], "01")
        assert_equal(variants.loc[7, "id"], "")
        assert_equal(variants.loc[7, "nalleles"], 7)
        assert_equal(variants.loc[7, "allele_ids"], "A,G,GT,GTT,GTTT,GTTTT,GTTTTT")
        assert_equal(variants.loc[7, "pos"], 8)
        assert_equal(variants.loc[7, "rsid"], "M8")

        n = variants.shape[0]
        assert_equal(variants.loc[n - 1, "chrom"], "01")
        assert_equal(variants.loc[n - 1, "id"], "")
        assert_equal(variants.loc[n - 1, "nalleles"], 2)
        assert_equal(variants.loc[n - 1, "allele_ids"], "A,G")
        assert_equal(variants.loc[n - 1, "pos"], 10)
        assert_equal(variants.loc[n - 1, "rsid"], "M10")

        assert_equal(samples.loc[0, "id"], "sample_0")
        assert_equal(samples.loc[3, "id"], "sample_3")

        G = bgen["genotype"][0, 0, :].compute()
        assert_allclose(G[:2], [1, 0])
        assert_(isnan(G[2]))

        G = bgen["genotype"][0, 1, :].compute()
        assert_allclose(G[:3], [1, 0, 0])

        G = bgen["genotype"][-1, -1, :].compute()
        assert_allclose(G[:5], [0, 0, 0, 1, 0])

        X = bgen["X"]

        assert_allclose(X[0].compute().sel(data="ploidy"), [1, 2, 2, 2])
        assert_allclose(X[-1].compute().sel(data="ploidy"), [4, 4, 4, 4])

        assert_allclose(
            X[:, 0].compute().sel(data="phased"), [0, 1, 1, 0, 1, 1, 1, 1, 0, 0]
        )

        X = X.compute()

        x = X.sel(sample=0, data="phased")
        assert_allclose(x.where(x == 1, drop=True).variant.values, [1, 2, 4, 5, 6, 7])


def test_bgen_reader_complex_sample_file():
    with example_files(["complex.23bits.bgen", "complex.sample"]) as filepaths:
        bgen = read_bgen(filepaths[0], sample_file=filepaths[1], verbose=False)
        variants = bgen["variants"]
        samples = bgen["samples"]
        assert_("genotype" in bgen)

        assert_equal(variants.loc[0, "chrom"], "01")
        assert_equal(variants.loc[0, "id"], "")
        assert_equal(variants.loc[0, "nalleles"], 2)
        assert_equal(variants.loc[0, "allele_ids"], "A,G")
        assert_equal(variants.loc[0, "pos"], 1)
        assert_equal(variants.loc[0, "rsid"], "V1")

        assert_equal(variants.loc[7, "chrom"], "01")
        assert_equal(variants.loc[7, "id"], "")
        assert_equal(variants.loc[7, "nalleles"], 7)
        assert_equal(variants.loc[7, "allele_ids"], "A,G,GT,GTT,GTTT,GTTTT,GTTTTT")
        assert_equal(variants.loc[7, "pos"], 8)
        assert_equal(variants.loc[7, "rsid"], "M8")

        n = variants.shape[0]
        assert_equal(variants.loc[n - 1, "chrom"], "01")
        assert_equal(variants.loc[n - 1, "id"], "")
        assert_equal(variants.loc[n - 1, "nalleles"], 2)
        assert_equal(variants.loc[n - 1, "allele_ids"], "A,G")
        assert_equal(variants.loc[n - 1, "pos"], 10)
        assert_equal(variants.loc[n - 1, "rsid"], "M10")

        assert_equal(samples.loc[0, "id"], "sample_0")
        assert_equal(samples.loc[3, "id"], "sample_3")
