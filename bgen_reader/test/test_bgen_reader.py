import os
import platform
from contextlib import contextmanager
from shutil import copyfile

import dask.dataframe as dd
import pytest
from dask.delayed import Delayed
from numpy import array, array_equal, isnan
from numpy.testing import assert_allclose, assert_equal
from pandas import Series

from bgen_reader import create_metafile, example_filepath, read_bgen


def test_bgen_samples_inside_bgen():
    data = read_bgen(example_filepath("haplotypes.bgen"), verbose=False)
    samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
    samples = Series(samples, dtype=str, name="id")
    assert all(data["samples"] == samples)


def test_bgen_samples_not_present():
    data = read_bgen(example_filepath("complex.23bits.no.samples.bgen"), verbose=False)
    samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
    samples = Series(samples, dtype=str, name="id")
    assert all(data["samples"] == samples)


def test_bgen_samples_specify_samples_file():
    data = read_bgen(
        example_filepath("complex.23bits.bgen"),
        samples_filepath=example_filepath("complex.sample"),
        verbose=False,
    )
    samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
    samples = Series(samples, dtype=str, name="id")
    assert all(data["samples"] == samples)


@pytest.mark.skipif(platform.system() != "Darwin", reason="only reliable on macos")
def test_bgen_samples_outside_bgen_unreadable(tmp_path):
    bgen_filepath = example_filepath("complex.23bits.bgen")
    samples_filepath = tmp_path / "complex.sample"
    copyfile(example_filepath("complex.sample"), samples_filepath)
    with noread_permission(samples_filepath):
        with pytest.raises(PermissionError):
            read_bgen(bgen_filepath, samples_filepath=samples_filepath, verbose=False)


@pytest.mark.skipif(platform.system() != "Darwin", reason="only reliable on macos")
def test_bgen_file_not_readable(tmp_path):
    filepath = tmp_path / "haplotypes.bgen"
    copyfile(example_filepath("haplotypes.bgen"), filepath)
    with noread_permission(filepath):
        with pytest.raises(PermissionError):
            read_bgen(filepath, verbose=False)


def test_bgen_file_dont_exist():
    with pytest.raises(FileNotFoundError):
        read_bgen("idontexist.bgen", verbose=False)


def test_metafile_not_provided():
    read_bgen(example_filepath("haplotypes.bgen"), verbose=False)


def test_metafile_provided_not_supported_anymore():
    with pytest.raises(RuntimeError):
        read_bgen(
            example_filepath("haplotypes.bgen"),
            metafile_filepath=example_filepath("haplotypes.bgen.metadata.valid"),
            verbose=False,
        )


def test_metafile_wrong_filepath():
    filepath = example_filepath("haplotypes.bgen")
    fp = "/omg/invalid/haplotypes.bgen.metafile_path"
    with pytest.raises(FileNotFoundError):
        with pytest.warns(UserWarning):
            read_bgen(filepath, metafile_filepath=fp, verbose=False)


@pytest.mark.skipif(platform.system() != "Darwin", reason="only reliable on macos")
def test_metafile_not_provided_no_permission_to_create(tmp_path):
    src = example_filepath("haplotypes.bgen")
    dst = tmp_path / "haplotypes.bgen"
    copyfile(src, dst)
    path = os.path.dirname(dst)
    with nowrite_permission(path):
        with pytest.warns(UserWarning):
            read_bgen(dst, verbose=False)


@contextmanager
def nowrite_permission(path):
    perm = os.stat(path).st_mode
    os.chmod(path, 0o555)
    try:
        yield
    finally:
        os.chmod(path, perm)


@contextmanager
def noread_permission(path):
    perm = os.stat(path).st_mode
    os.chmod(path, 0o333)
    try:
        yield
    finally:
        os.chmod(path, perm)


def test_bgen_reader_lazy_types():
    bgen = read_bgen(example_filepath("haplotypes.bgen"), verbose=False)
    assert isinstance(bgen["genotype"][0], Delayed)
    assert isinstance(bgen["variants"], dd.DataFrame)


def test_bgen_reader_phased_genotype():
    filepath = example_filepath("haplotypes.bgen")
    bgen = read_bgen(filepath, verbose=False)
    variants = bgen["variants"]
    samples = bgen["samples"]

    v = variants.loc[0].compute()
    assert_equal(v["chrom"].values[0], "1")
    assert_equal(v["id"].values[0], "SNP1")
    assert_equal(v["nalleles"].values[0], 2)
    assert_equal(v["allele_ids"].values[0], "A,G")
    assert_equal(v["pos"].values[0], 1)
    assert_equal(v["rsid"].values[0], "RS1")

    v = variants.loc[2].compute()
    assert_equal(v["chrom"].values[0], "1")
    assert_equal(v["id"].values[0], "SNP3")
    assert_equal(v["nalleles"].values[0], 2)
    assert_equal(v["allele_ids"].values[0], "A,G")
    assert_equal(v["pos"].values[0], 3)
    assert_equal(v["rsid"].values[0], "RS3")

    assert_equal(samples.loc[0], "sample_0")
    assert_equal(samples.loc[2], "sample_2")

    n = samples.shape[0]
    assert_equal(samples.loc[n - 1], "sample_3")

    g = bgen["genotype"][0].compute()
    assert_allclose(g["probs"][0], [1.0, 0.0, 1.0, 0.0])
    k = len(variants)
    n = len(samples)
    g = bgen["genotype"][k - 1].compute()
    assert_allclose(g["probs"][n - 1], [1.0, 0.0, 0.0, 1.0])


def test_bgen_reader_variants_info():
    filepath = example_filepath("example.32bits.bgen")
    bgen = read_bgen(filepath, verbose=False)
    variants = bgen["variants"]
    samples = bgen["samples"]
    assert "genotype" in bgen

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

    g = bgen["genotype"][0].compute()["probs"]
    assert all(isnan(g[0, :]))

    g = bgen["genotype"][0].compute()["probs"]
    a = [0.027802362811705648, 0.00863673794284387, 0.9635608992454505]
    assert_allclose(g[1, :], a)

    b = [
        0.97970582847010945215516,
        0.01947019668749305418287,
        0.00082397484239749366197,
    ]
    g = bgen["genotype"][1].compute()["probs"]
    assert_allclose(g[2, :], b)


def test_bgen_reader_without_metadata():
    filepath = example_filepath("example.32bits.bgen")
    bgen = read_bgen(filepath, verbose=False)
    variants = bgen["variants"].compute()
    samples = bgen["samples"]
    assert "genotype" in bgen
    assert_equal(variants.loc[7, "allele_ids"], "A,G")
    n = samples.shape[0]
    assert_equal(samples.loc[n - 1], "sample_500")


def test_bgen_reader_with_wrong_metadata_file():
    filepath = example_filepath("example.32bits.bgen")
    filepath.touch()
    metafile_filepath = example_filepath("wrong.metadata")
    metafile_filepath.touch()  # make sure that the metafile has a later timestamp (otherwise, it might be re-created)
    with pytest.raises(RuntimeError):
        read_bgen(filepath, verbose=False, metafile_filepath=metafile_filepath)


def test_bgen_reader_with_nonexistent_metadata_file():
    filepath = example_filepath("example.32bits.bgen")
    folder = os.path.dirname(filepath)
    metafile_filepath = os.path.join(folder, "nonexistent.metadata")

    with pytest.raises(FileNotFoundError):
        with pytest.warns(UserWarning):
            read_bgen(filepath, verbose=False, metafile_filepath=metafile_filepath)


def test_bgen_reader_file_notfound():
    with pytest.raises(FileNotFoundError):
        read_bgen("/1/2/3/example.32bits.bgen", verbose=False)


def test_create_metadata_file(tmp_path):
    filepath = example_filepath("example.32bits.bgen")
    metafile_filepath = tmp_path / (filepath.name + ".metadata")

    try:
        create_metafile(filepath, metafile_filepath, verbose=False)
        assert os.path.exists(metafile_filepath)
    finally:
        if os.path.exists(metafile_filepath):
            os.remove(metafile_filepath)


def test_bgen_reader_complex():
    filepath = example_filepath("complex.23bits.bgen")
    bgen = read_bgen(filepath, verbose=False)
    variants = bgen["variants"].compute()
    samples = bgen["samples"]
    assert "genotype" in bgen

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

    g = bgen["genotype"][0].compute()["probs"][0]
    assert_allclose(g[:2], [1, 0])
    assert isnan(g[2])

    g = bgen["genotype"][0].compute()["probs"][1]
    assert_allclose(g[:3], [1, 0, 0])

    g = bgen["genotype"][-1].compute()["probs"][-1]
    assert_allclose(g[:5], [0, 0, 0, 1, 0])

    ploidy = bgen["genotype"][0].compute()["ploidy"]
    assert_allclose(ploidy, [1, 2, 2, 2])
    ploidy = bgen["genotype"][-1].compute()["ploidy"]
    assert_allclose(ploidy, [4, 4, 4, 4])

    nvariants = len(variants)
    phased = [bgen["genotype"][i].compute()["phased"] for i in range(nvariants)]
    phased = array(phased)
    assert_equal(phased.dtype.name, "bool")
    ideal = array([False, True, True, False, True, True, True, True, False, False])
    assert array_equal(phased, ideal)


def test_bgen_reader_complex_sample_file():
    bgen = read_bgen(
        example_filepath("complex.23bits.bgen"),
        samples_filepath=example_filepath("complex.sample"),
        verbose=False,
    )
    variants = bgen["variants"].compute()
    samples = bgen["samples"]
    assert "genotype" in bgen

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

    ploidy = bgen["genotype"][2].compute()["ploidy"]
    missing = bgen["genotype"][2].compute()["missing"]
    nvariants = len(variants)
    phased = [bgen["genotype"][i].compute()["phased"] for i in range(nvariants)]
    assert_allclose(ploidy, [1, 2, 2, 2])
    assert_allclose(missing, [0, 0, 0, 0])
    assert_allclose(phased, [0, 1, 1, 0, 1, 1, 1, 1, 0, 0])


if __name__ == "__main__":
    pytest.main([__file__])
