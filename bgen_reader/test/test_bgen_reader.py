from __future__ import unicode_literals

import os
import stat
import sys
from contextlib import contextmanager

import pytest
from numpy import isnan
from numpy.testing import assert_, assert_allclose, assert_equal
from pandas import Series

from bgen_reader import create_metafile, example_files, read_bgen

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def test_bgen_samples_inside_bgen():
    with example_files("haplotypes.bgen") as filepath:
        data = read_bgen(filepath, verbose=False)
        samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
        samples = Series(samples, dtype=str, name="id")
        assert_(all(data["samples"] == samples))


def test_bgen_samples_not_present():
    with example_files("complex.23bits.no.samples.bgen") as filepath:
        data = read_bgen(filepath, verbose=False)
        samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
        samples = Series(samples, dtype=str, name="id")
        assert_(all(data["samples"] == samples))


def test_bgen_samples_specify_samples_file():
    with example_files(["complex.23bits.bgen", "complex.sample"]) as filepaths:
        data = read_bgen(filepaths[0], samples_filepath=filepaths[1], verbose=False)
        samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
        samples = Series(samples, dtype=str, name="id")
        assert_(all(data["samples"] == samples))


def test_bgen_samples_outside_bgen_unreadable():
    with example_files(["complex.23bits.bgen", "complex.sample"]) as filepaths:
        with noread_permission(filepaths[1]):
            with pytest.raises(PermissionError):
                read_bgen(filepaths[0], samples_filepath=filepaths[1], verbose=False)


def test_bgen_file_not_readable():
    with example_files("haplotypes.bgen") as filepath:
        with pytest.raises(PermissionError):
            with noread_permission(filepath):
                read_bgen(filepath, verbose=False)


def test_bgen_file_dont_exist():
    with pytest.raises(FileNotFoundError):
        read_bgen("idontexist.bgen", verbose=False)


def test_metafile_not_provided():
    with example_files("haplotypes.bgen") as filepath:
        read_bgen(filepath, verbose=False)


def test_metafile_provided():
    filenames = ["haplotypes.bgen", "haplotypes.bgen.metadata.valid"]
    with example_files(filenames) as filepaths:
        read_bgen(filepaths[0], metafile_filepath=filepaths[1], verbose=False)


def test_metafile_wrong_filepath():
    with example_files("haplotypes.bgen") as filepath:
        fp = "/omg/invalid/haplotypes.bgen.metafile_path"
        with pytest.raises(FileNotFoundError):
            with pytest.warns(UserWarning):
                read_bgen(filepath, metafile_filepath=fp, verbose=False)


def test_metafile_not_provided_no_permission_to_create():
    with example_files("haplotypes.bgen") as filepath:
        path = os.path.dirname(filepath)
        with nowrite_permission(path):
            with pytest.warns(UserWarning):
                read_bgen(filepath, verbose=False)


@contextmanager
def nowrite_permission(path):
    perm = os.stat(path).st_mode
    os.chmod(path, stat.S_IXUSR | stat.S_IXGRP | stat.S_IRUSR | stat.S_IRGRP)
    try:
        yield
    finally:
        os.chmod(path, perm)


@contextmanager
def noread_permission(path):
    perm = os.stat(path).st_mode
    os.chmod(path, stat.S_IXUSR | stat.S_IXGRP)
    try:
        yield
    finally:
        os.chmod(path, perm)


def test_bgen_reader_phased_genotype():
    with example_files("haplotypes.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"]
        samples = bgen["samples"]

        v = variants.loc[0].compute()
        assert_equal(v["chrom"].item(), "1")
        assert_equal(v["id"].item(), "SNP1")
        assert_equal(v["nalleles"].item(), 2)
        assert_equal(v["allele_ids"].item(), "A,G")
        assert_equal(v["pos"].item(), 1)
        assert_equal(v["rsid"].item(), "RS1")

        v = variants.loc[2].compute()
        assert_equal(v["chrom"].item(), "1")
        assert_equal(v["id"].item(), "SNP3")
        assert_equal(v["nalleles"].item(), 2)
        assert_equal(v["allele_ids"].item(), "A,G")
        assert_equal(v["pos"].item(), 3)
        assert_equal(v["rsid"].item(), "RS3")

        assert_equal(samples.loc[0], "sample_0")
        assert_equal(samples.loc[2], "sample_2")

        n = samples.shape[0]
        assert_equal(samples.loc[n - 1], "sample_3")

        assert_allclose(bgen["genotype"][0]["probs"][0], [1.0, 0.0, 1.0, 0.0])
        k = len(variants)
        n = len(samples)
        assert_allclose(bgen["genotype"][k - 1]["probs"][n - 1], [1.0, 0.0, 0.0, 1.0])


def test_bgen_reader_variants_info():
    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"]
        samples = bgen["samples"]
        assert_("genotype" in bgen)

        variants = variants.compute()
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

        assert_equal(samples.loc[0], "sample_001")
        assert_equal(samples.loc[7], "sample_008")

        n = samples.shape[0]
        assert_equal(samples.loc[n - 1], "sample_500")

        g = bgen["genotype"][0]["probs"]
        assert_(all(isnan(g[0, :])))

        g = bgen["genotype"][0]["probs"]
        a = [0.027802362811705648, 0.00863673794284387, 0.9635608992454505]
        assert_allclose(g[1, :], a)

        b = [
            0.97970582847010945215516,
            0.01947019668749305418287,
            0.00082397484239749366197,
        ]
        g = bgen["genotype"][1]["probs"]
        assert_allclose(g[2, :], b)


def _test_bgen_reader_phased_genotype():
    with example_files("haplotypes.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"].compute()
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

        assert_equal(samples.loc[0], "sample_0")
        assert_equal(samples.loc[2], "sample_2")

        n = samples.shape[0]
        assert_equal(samples.loc[n - 1], "sample_3")

        g = bgen["genotype"][0]["probs"]
        a = [1.0, 0.0, 1.0, 0.0]
        assert_allclose(g[0, :], a)

        k = len(variants)
        n = len(samples)
        a = [1.0, 0.0, 0.0, 1.0]
        g = bgen["genotype"][k - 1]["probs"]
        assert_allclose(g[n - 1, :], a)


def test_bgen_reader_without_metadata():
    with example_files("example.32bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"].compute()
        samples = bgen["samples"]
        assert_("genotype" in bgen)
        assert_equal(variants.loc[7, "allele_ids"], "A,G")
        n = samples.shape[0]
        assert_equal(samples.loc[n - 1], "sample_500")


def test_bgen_reader_with_wrong_metadata_file():
    with example_files(["example.32bits.bgen", "wrong.metadata"]) as filepaths:
        with pytest.raises(RuntimeError):
            read_bgen(filepaths[0], verbose=False, metafile_filepath=filepaths[1])


def test_bgen_reader_with_nonexistent_metadata_file():
    with example_files("example.32bits.bgen") as filepath:
        folder = os.path.dirname(filepath)
        metafile_filepath = os.path.join(folder, "nonexistent.metadata")

        with pytest.raises(FileNotFoundError):
            read_bgen(filepath, verbose=False, metafile_filepath=metafile_filepath)


def test_bgen_reader_file_notfound():
    with pytest.raises(FileNotFoundError):
        read_bgen("/1/2/3/example.32bits.bgen", verbose=False)


# def test_create_metadata_file():
#     with example_files("example.32bits.bgen") as filepath:
#         folder = os.path.dirname(filepath)
#         metafile_filepath = os.path.join(folder, filepath + ".metadata")

#         create_metafile(filepath, metafile_filepath, verbose=False)
#         assert_(os.path.exists(metafile_filepath))
#         os.remove(metafile_filepath)

# @pytest.mark.skipif(sys.version_info < (3, 6), reason="requires python3.6 or higher")
# def test_create_metadata_file_bytes():
#     with example_files("example.32bits.bgen") as filepath:
#         folder = os.path.dirname(filepath)
#         metafile_filepath = os.path.join(folder, filepath + ".metadata")

#         create_metafile(filepath, metafile_filepath.encode(), verbose=False)
#         assert_(os.path.exists(metafile_filepath))
#         os.remove(metafile_filepath)


def test_bgen_reader_complex():
    with example_files("complex.23bits.bgen") as filepath:
        bgen = read_bgen(filepath, verbose=False)
        variants = bgen["variants"].compute()
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

        assert_equal(samples.loc[0], "sample_0")
        assert_equal(samples.loc[3], "sample_3")

        g = bgen["genotype"][0]["probs"][0]
        assert_allclose(g[:2], [1, 0])
        assert_(isnan(g[2]))

        g = bgen["genotype"][0]["probs"][1]
        assert_allclose(g[:3], [1, 0, 0])

        g = bgen["genotype"][-1]["probs"][-1]
        assert_allclose(g[:5], [0, 0, 0, 1, 0])

#         X = bgen["X"]

#         assert_allclose(X[0].compute().sel(data="ploidy"), [1, 2, 2, 2])
#         assert_allclose(X[-1].compute().sel(data="ploidy"), [4, 4, 4, 4])

#         assert_allclose(
#             X[:, 0].compute().sel(data="phased"), [0, 1, 1, 0, 1, 1, 1, 1, 0, 0]
#         )

#         X = X.compute()

#         x = X.sel(sample=0, data="phased")
#         assert_allclose(x.where(x == 1, drop=True).variant.values, [1, 2, 4, 5, 6, 7])

# def test_bgen_reader_complex2():
#     with example_files("complex.23bits.bgen") as filepath:
#         bgen = read_bgen(filepath, verbose=False)
#         variants = bgen["variants"]
#         samples = bgen["samples"]

#         assert_equal(variants.loc[0, "chrom"].compute().item(), "01")
#         assert_equal(variants.loc[0, "id"].compute().item(), "")
#         assert_equal(variants.loc[0, "nalleles"].compute().item(), 2)
#         assert_equal(variants.loc[0, "allele_ids"].compute().item(), "A,G")
#         assert_equal(variants.loc[0, "pos"].compute().item(), 1)
#         assert_equal(variants.loc[0, "rsid"].compute().item(), "V1")

#         assert_equal(variants.loc[7, "chrom"].compute().item(), "01")
#         assert_equal(variants.loc[7, "id"].compute().item(), "")
#         assert_equal(variants.loc[7, "nalleles"].compute().item(), 7)
#         assert_equal(
#             variants.loc[7, "allele_ids"].compute().item(),
#             "A,G,GT,GTT,GTTT,GTTTT,GTTTTT",
#         )
#         assert_equal(variants.loc[7, "pos"].compute().item(), 8)
#         assert_equal(variants.loc[7, "rsid"].compute().item(), "M8")

#         n = variants.shape[0].compute()
#         assert_equal(variants.loc[n - 1, "chrom"].compute().item(), "01")
#         assert_equal(variants.loc[n - 1, "id"].compute().item(), "")
#         assert_equal(variants.loc[n - 1, "nalleles"].compute().item(), 2)
#         assert_equal(variants.loc[n - 1, "allele_ids"].compute().item(), "A,G")
#         assert_equal(variants.loc[n - 1, "pos"].compute().item(), 10)
#         assert_equal(variants.loc[n - 1, "rsid"].compute().item(), "M10")

#         assert_equal(samples.loc[0], "sample_0")
#         assert_equal(samples.loc[3], "sample_3")

#         G = bgen["variants"]["genotype"].loc[0].compute().item().compute()
#         assert_allclose(G[0, :2], [1, 0])
#         assert_(isnan(G[0, 2]))

#         G = bgen["variants"]["genotype"].loc[0].compute().item().compute()
#         assert_allclose(G[2, :3], [1, 0, 0])

#         G = bgen["variants"]["genotype"].compute().values[-1].compute()
#         assert_allclose(G[-1, :5], [0, 0, 0, 1, 0])

#         assert_allclose(
#             bgen["variants"]["genotype"].loc[0].compute().item().compute().ploidy,
#             [1, 2, 2, 2],
#         )
#         assert_allclose(
#             bgen["variants"]["genotype"].compute().values[-1].compute().ploidy,
#             [4, 4, 4, 4],
#         )

#         assert_allclose(
#             [bgen["variants"]["genotype"].loc[i].compute().item().compute().phased for i in range(10)],
#             [0, 1, 1, 0, 1, 1, 1, 1, 0, 0],
#         )

#         # TODO: fix it
#         # x = X.sel(sample=0, data="phased")
#         # assert_allclose(x.where(x == 1, drop=True).variant.values, [1, 2, 4, 5, 6, 7])

# def test_bgen_reader_complex_sample_file():
#     with example_files(["complex.23bits.bgen", "complex.sample"]) as filepaths:
#         bgen = read_bgen(filepaths[0], sample_file=filepaths[1], verbose=False)
#         variants = bgen["variants"]
#         samples = bgen["samples"]
#         assert_("genotype" in bgen)

#         assert_equal(variants.loc[0, "chrom"], "01")
#         assert_equal(variants.loc[0, "id"], "")
#         assert_equal(variants.loc[0, "nalleles"], 2)
#         assert_equal(variants.loc[0, "allele_ids"], "A,G")
#         assert_equal(variants.loc[0, "pos"], 1)
#         assert_equal(variants.loc[0, "rsid"], "V1")

#         assert_equal(variants.loc[7, "chrom"], "01")
#         assert_equal(variants.loc[7, "id"], "")
#         assert_equal(variants.loc[7, "nalleles"], 7)
#         assert_equal(variants.loc[7, "allele_ids"], "A,G,GT,GTT,GTTT,GTTTT,GTTTTT")
#         assert_equal(variants.loc[7, "pos"], 8)
#         assert_equal(variants.loc[7, "rsid"], "M8")

#         n = variants.shape[0]
#         assert_equal(variants.loc[n - 1, "chrom"], "01")
#         assert_equal(variants.loc[n - 1, "id"], "")
#         assert_equal(variants.loc[n - 1, "nalleles"], 2)
#         assert_equal(variants.loc[n - 1, "allele_ids"], "A,G")
#         assert_equal(variants.loc[n - 1, "pos"], 10)
#         assert_equal(variants.loc[n - 1, "rsid"], "M10")

#         assert_equal(samples.loc[0, "id"], "sample_0")
#         assert_equal(samples.loc[3, "id"], "sample_3")

# def test_bgen_reader_complex_sample_file2():
#     with example_files(["complex.23bits.bgen", "complex.sample"]) as filepaths:
#         bgen = read_bgen(filepaths[0], sample_file=filepaths[1], verbose=False)
#         variants = bgen["variants"]
#         samples = bgen["samples"]

#         assert_equal(variants.loc[0, "chrom"].compute().item(), "01")
#         assert_equal(variants.loc[0, "id"].compute().item(), "")
#         assert_equal(variants.loc[0, "nalleles"].compute().item(), 2)
#         assert_equal(variants.loc[0, "allele_ids"].compute().item(), "A,G")
#         assert_equal(variants.loc[0, "pos"].compute().item(), 1)
#         assert_equal(variants.loc[0, "rsid"].compute().item(), "V1")

#         assert_equal(variants.loc[7, "chrom"].compute().item(), "01")
#         assert_equal(variants.loc[7, "id"].compute().item(), "")
#         assert_equal(variants.loc[7, "nalleles"].compute().item(), 7)
#         assert_equal(variants.loc[7, "allele_ids"].compute().item(), "A,G,GT,GTT,GTTT,GTTTT,GTTTTT")
#         assert_equal(variants.loc[7, "pos"].compute().item(), 8)
#         assert_equal(variants.loc[7, "rsid"].compute().item(), "M8")

#         n = variants.shape[0].compute()
#         assert_equal(variants["chrom"].compute().values[n - 1], "01")
#         assert_equal(variants["id"].compute().values[n - 1], "")
#         assert_equal(variants["nalleles"].compute().values[n - 1], 2)
#         assert_equal(variants["allele_ids"].compute().values[n - 1], "A,G")
#         assert_equal(variants["pos"].compute().values[n - 1], 10)
#         assert_equal(variants["rsid"].compute().values[n - 1], "M10")

#         assert_equal(samples.loc[0], "sample_0")
#         assert_equal(samples.loc[3], "sample_3")

# def test_bgen_reader_too_small_chunk_size():
#     with example_files("complex.23bits.bgen") as filepath:
#         read_bgen(filepath, size=1e-10, verbose=False)

# def test_bgen_reader_too_small_chunk_size2():
#     with example_files("complex.23bits.bgen") as filepath:
#         read_bgen(filepath, size=1e-10, verbose=False)

# def test_bgen_reader_concat_chunks():
#     with example_files("complex.23bits.bgen") as filepath:
#         bgen = read_bgen(filepath, size=0.0001, verbose=False)

#         # G = bgen["genotype"][0, 0, :].compute()
#         # print(type(G))

#         # G = bgen["genotype"][0, [0, 1]].compute()
#         # print(type(G))

#         # G = bgen["genotype"][[0, 1]].compute()
#         # print(type(G))

#         # G = bgen["genotype"][[0]].compute()
#         # print(type(G))

#         # G = bgen["genotype"].compute()
#         # print(type(G))

#         # print(bgen["genotype"].chunks)
#         # assert_allclose(G[:2], [1, 0])
#         # assert_(isnan(G[2]))

#         # G = bgen["genotype"][0, 1, :].compute()
#         # assert_allclose(G[:3], [1, 0, 0])

#         # G = bgen["genotype"][-1, -1, :].compute()
#         # assert_allclose(G[:5], [0, 0, 0, 1, 0])

#         # X = bgen["X"]

#         # assert_allclose(X[0].compute().sel(data="ploidy"), [1, 2, 2, 2])
#         # assert_allclose(X[-1].compute().sel(data="ploidy"), [4, 4, 4, 4])

#         # assert_allclose(
#         #     X[:, 0].compute().sel(data="phased"), [0, 1, 1, 0, 1, 1, 1, 1, 0, 0]
#         # )

#         # X = X.compute()

#         # x = X.sel(sample=0, data="phased")
#         # assert_allclose(x.where(x == 1, drop=True).variant.values, [1, 2, 4, 5, 6, 7])
