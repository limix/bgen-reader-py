import os
from contextlib import contextmanager
from shutil import copyfile
import platform
from pathlib import Path
import numpy as np

import pytest
from numpy import array, array_equal, isnan
from numpy.testing import assert_allclose, assert_equal
from pandas import Series

from bgen_reader import open_bgen #!!!cmk switch to relative
from bgen_reader import example_filepath
from bgen_reader._environment import BGEN_READER_CACHE_HOME
from bgen_reader.test.write_random import _write_random

def example_filepath2(filename):
    filepath = example_filepath(filename)
    metadata2_path = open_bgen._metadatapath_from_filename(filepath)
    if metadata2_path.exists():
        metadata2_path.unlink()
    return filepath

def test_bgen_samples_inside_bgen():
    data = open_bgen(example_filepath2("haplotypes.bgen"), verbose=False)
    samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
    assert all(data.samples == samples)


def test_bgen_samples_not_present():
    data = open_bgen(example_filepath2("complex.23bits.no.samples.bgen"), verbose=False)
    samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
    assert all(data.samples == samples)


def test_bgen_samples_specify_samples_file():
    data = open_bgen(
        example_filepath2("complex.23bits.bgen"),
        samples_filepath=example_filepath("complex.sample"),
        verbose=False,
    )
    samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
    assert all(data.samples == samples)


@pytest.mark.skipif(platform.system() != "Darwin", reason="only reliable on macos")
def test_bgen_samples_outside_bgen_unreadable(tmp_path):
    bgen_filepath = example_filepath2("complex.23bits.bgen")
    samples_filepath = tmp_path / "complex.sample"
    copyfile(example_filepath("complex.sample"), samples_filepath)
    with noread_permission(samples_filepath):
        with pytest.raises(PermissionError):
            open_bgen(bgen_filepath, samples_filepath=samples_filepath, verbose=False)


@pytest.mark.skipif(platform.system() != "Darwin", reason="only reliable on macos")
def test_bgen_file_not_readable(tmp_path):
    filepath = tmp_path / "haplotypes.bgen"
    copyfile(example_filepath2("haplotypes.bgen"), filepath)
    with noread_permission(filepath):
        with pytest.raises(PermissionError):
            open_bgen(filepath, verbose=False)


def test_bgen_file_dont_exist():
    with pytest.raises(FileNotFoundError):
        open_bgen("idontexist.bgen", verbose=False)


def test_metafile_not_provided():
    open_bgen(example_filepath2("haplotypes.bgen"), verbose=False)


@pytest.mark.skipif(platform.system() != "Darwin", reason="only reliable on macos")
def test_metafile_not_provided_no_permission_to_create(tmp_path):
    src = example_filepath2("haplotypes.bgen")
    dst = tmp_path / "haplotypes.bgen"
    copyfile(src, dst)
    path = os.path.dirname(dst)
    with nowrite_permission(path):
        with pytest.warns(UserWarning):
            open_bgen(dst, verbose=False)


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

def test_open_bgen_phased_genotype():
    filepath = example_filepath2("haplotypes.bgen")
    bgen2 = open_bgen(filepath, verbose=False)

    assert_equal(bgen2.chromosomes[0], "1")
    assert_equal(bgen2.ids[0], "SNP1")
    assert_equal(bgen2.nalleles[0], 2)
    assert_equal(bgen2.allele_ids[0], "A,G")
    assert_equal(bgen2.positions[0], 1)
    assert_equal(bgen2.rsids[0], "RS1")

    assert_equal(bgen2.chromosomes[2], "1")
    assert_equal(bgen2.ids[2], "SNP3")
    assert_equal(bgen2.nalleles[2], 2)
    assert_equal(bgen2.allele_ids[2], "A,G")
    assert_equal(bgen2.positions[2], 3)
    assert_equal(bgen2.rsids[2], "RS3")

    assert_equal(bgen2.samples[0], "sample_0")
    assert_equal(bgen2.samples[2], "sample_2")
    assert_equal(bgen2.samples[-1], "sample_3")

    g = bgen2.read((0,0))
    assert_allclose(g[0,0,:], [1.0, 0.0, 1.0, 0.0])
    g = bgen2.read((-1,-1))
    assert_allclose(g[0,0,:], [1.0, 0.0, 0.0, 1.0])


def test_open_bgen_variants_info():
    filepath = example_filepath2("example.32bits.bgen")
    bgen2 = open_bgen(filepath, verbose=False)

    assert_equal(bgen2.chromosomes[0], "01")
    assert_equal(bgen2.ids[0], "SNPID_2")
    assert_equal(bgen2.nalleles[0], 2)
    assert_equal(bgen2.allele_ids[0], "A,G")
    assert_equal(bgen2.positions[0], 2000)
    assert_equal(bgen2.rsids[0], "RSID_2")

    assert_equal(bgen2.chromosomes[7], "01")
    assert_equal(bgen2.ids[7], "SNPID_9")
    assert_equal(bgen2.nalleles[7], 2)
    assert_equal(bgen2.allele_ids[7], "A,G")
    assert_equal(bgen2.positions[7], 9000)
    assert_equal(bgen2.rsids[7], "RSID_9")

    assert_equal(bgen2.chromosomes[-1], "01")
    assert_equal(bgen2.ids[-1], "SNPID_200")
    assert_equal(bgen2.nalleles[-1], 2)
    assert_equal(bgen2.allele_ids[-1], "A,G")
    assert_equal(bgen2.positions[-1], 100001)
    assert_equal(bgen2.rsids[-1], "RSID_200")

    assert_equal(bgen2.samples[0], "sample_001")
    assert_equal(bgen2.samples[7], "sample_008")
    assert_equal(bgen2.samples[-1], "sample_500")

    g = bgen2.read((0,0))
    assert all(isnan(g[0,0,:]))

    g = bgen2.read((1,0))
    a = [0.027802362811705648, 0.00863673794284387, 0.9635608992454505]
    assert_allclose(g[0,0,:], a)

    b = [
        0.97970582847010945215516,
        0.01947019668749305418287,
        0.00082397484239749366197,
    ]
    g = bgen2.read((2,1))
    assert_allclose(g[0,0,:], b)

def test_to_improve_coverage():
    filepath = example_filepath2("example.32bits.bgen")
    bgen2 = open_bgen(filepath, verbose=False) #Creates metadata2.npz file
    assert_equal(bgen2.ncombinations[-1],3)
    assert_equal(bgen2.phased[-1], False)
    with open_bgen(filepath) as bgen2: #Reuses metadata2.npz file
        assert_equal(str(bgen2),"open_bgen('{0}')".format(filepath.name))
        assert_equal(bgen2.nsamples, 500)
        assert_equal(bgen2.nvariants, 199)
        assert_equal(bgen2.shape, (500,199,3))
        assert_equal(bgen2.ids[-1], "SNPID_200")
        assert_equal(bgen2.rsids[-1], "RSID_200")
        assert_equal(bgen2.chromosomes[-1], "01")
        assert_equal(bgen2.positions[-1], 100001)
        assert_equal(bgen2.nalleles[-1], 2)
        assert_equal(bgen2.allele_ids[-1], "A,G")
        assert_equal(bgen2.ncombinations[-1],3)
        assert_equal(bgen2.phased[-1], False)
        assert_equal(bgen2.samples[-1], "sample_500")

        b = [
            0.97970582847010945215516,
            0.01947019668749305418287,
            0.00082397484239749366197,
        ]
        g = bgen2.read((2,1))
        assert_allclose(g[0,0,:], b)

        g = bgen2.read()
        assert_allclose(g[2,1,:], b)

def test_bigfile():
    nsamples = 2500
    nvariants = 500*1000
    bits=16
    test_data_folder = BGEN_READER_CACHE_HOME / "test_data"
    filepath = test_data_folder / "{0}x{1}.{2}bits.bgen".format(nsamples,nvariants,bits)
    if not filepath.exists():
        _write_random(filepath,nsamples,nvariants,bits=bits,cleanup_temp_files=False) #!!!cmk change to cleanup_temp_files=True (the default)
    #!!!cmk remove metadata2.npz file
    with open_bgen(filepath,verbose=True) as bgen2:
        assert bgen2.nsamples == nsamples #!!!cmk use other asserts?
        assert bgen2.nvariants == nvariants
        val = bgen2.read(-1)
        assert val.shape == (nsamples,1,3)
        mean = np.nanmean(val)
        assert mean!=mean or (0 <= mean and mean <= 1)


#!!!cmk why is there also a (no '_') test_open_bgen_phased_genotype?
def _test_open_bgen_phased_genotype():
    filepath = example_filepath2("haplotypes.bgen")
    remove_any_metadata2_file(filepath)
    bgen2 = open_bgen(filepath, verbose=False)

    assert_equal(bgen2.chromosomes[0], "1")
    assert_equal(bgen2.ids[0], "SNP1")
    assert_equal(bgen2.nalleles[0], 2)
    assert_equal(bgen2.allele_ids[0], "A,G")
    assert_equal(bgen2.positions[0], 1)
    assert_equal(bgen2.rsids[0], "RS1")

    assert_equal(bgen2.chromosomes[2], "1")
    assert_equal(bgen2.ids[2], "SNP3")
    assert_equal(bgen2.nalleles[2], 2)
    assert_equal(bgen2.allele_ids[2], "A,G")
    assert_equal(bgen2.positions[2], 3)
    assert_equal(bgen2.rsids[2], "RS3")

    assert_equal(bgen2.samples[0], "sample_0")
    assert_equal(bgen2.samples[2], "sample_2")
    assert_equal(bgen2.samples[-1], "sample_3")

    g = bgen2.read((0,0))
    a = [1.0, 0.0, 1.0, 0.0]
    assert_allclose(g[0,0,:], a)

    a = [1.0, 0.0, 0.0, 1.0]
    g = bgen2.read((-1,-1))
    assert_allclose(g[0,0,:], a)


def test_open_bgen_without_metadata():
    filepath = example_filepath2("example.32bits.bgen")
    bgen2 = open_bgen(filepath, verbose=False)
    assert_equal(bgen2.allele_ids[7], "A,G")
    assert_equal(bgen2.samples[-1], "sample_500")

def test_open_bgen_file_notfound():
    with pytest.raises(FileNotFoundError):
        open_bgen("/1/2/3/example.32bits.bgen", verbose=False)


def test_open_bgen_complex():
    filepath = example_filepath2("complex.23bits.bgen")
    bgen2 = open_bgen(filepath, verbose=False)

    assert_equal(bgen2.chromosomes[0], "01")
    assert_equal(bgen2.ids[0], "")
    assert_equal(bgen2.nalleles[0], 2)
    assert_equal(bgen2.allele_ids[0], "A,G")
    assert_equal(bgen2.positions[0], 1)
    assert_equal(bgen2.rsids[0], "V1")

    assert_equal(bgen2.chromosomes[7], "01")
    assert_equal(bgen2.ids[7], "")
    assert_equal(bgen2.nalleles[7], 7)
    assert_equal(bgen2.allele_ids[7], "A,G,GT,GTT,GTTT,GTTTT,GTTTTT")
    assert_equal(bgen2.positions[7], 8)
    assert_equal(bgen2.rsids[7], "M8")

    assert_equal(bgen2.chromosomes[-1], "01")
    assert_equal(bgen2.ids[-1], "")
    assert_equal(bgen2.nalleles[-1], 2)
    assert_equal(bgen2.allele_ids[-1], "A,G")
    assert_equal(bgen2.positions[-1], 10)
    assert_equal(bgen2.rsids[-1], "M10")

    assert_equal(bgen2.samples[0], "sample_0")
    assert_equal(bgen2.samples[3], "sample_3")

    g = bgen2.read((0,0))
    assert_allclose(g[0,0,:2], [1, 0])
    assert isnan(g[0,0,2])

    g = bgen2.read((1,0))
    assert_allclose(g[0,0,:3], [1, 0, 0])

    g = bgen2.read((-1,-1))
    assert_allclose(g[0,0,:5], [0, 0, 0, 1, 0])

    ploidy = bgen2.read(0,return_probabilities=False,return_ploidies=True)
    assert_allclose(ploidy[:,0], [1, 2, 2, 2])
    ploidy = bgen2.read(-1,return_probabilities=False,return_ploidies=True)
    assert_allclose(ploidy[:,0], [4, 4, 4, 4])

    assert_equal(bgen2.phased.dtype.name, "bool")
    ideal = array([False, True, True, False, True, True, True, True, False, False])
    assert array_equal(bgen2.phased, ideal)


def test_open_bgen_complex_sample_file():
    bgen2 = open_bgen(
        example_filepath2("complex.23bits.bgen"),
        samples_filepath=example_filepath("complex.sample"),
        verbose=False,
    )

    assert_equal(bgen2.chromosomes[0], "01")
    assert_equal(bgen2.ids[0], "")
    assert_equal(bgen2.nalleles[0], 2)
    assert_equal(bgen2.allele_ids[0], "A,G")
    assert_equal(bgen2.positions[0], 1)
    assert_equal(bgen2.rsids[0], "V1")

    assert_equal(bgen2.chromosomes[7], "01")
    assert_equal(bgen2.ids[7], "")
    assert_equal(bgen2.nalleles[7], 7)
    assert_equal(bgen2.allele_ids[7], "A,G,GT,GTT,GTTT,GTTTT,GTTTTT")
    assert_equal(bgen2.positions[7], 8)
    assert_equal(bgen2.rsids[7], "M8")

    assert_equal(bgen2.chromosomes[-1], "01")
    assert_equal(bgen2.ids[-1], "")
    assert_equal(bgen2.nalleles[-1], 2)
    assert_equal(bgen2.allele_ids[-1], "A,G")
    assert_equal(bgen2.positions[-1], 10)
    assert_equal(bgen2.rsids[-1], "M10")

    assert_equal(bgen2.samples[0], "sample_0")
    assert_equal(bgen2.samples[3], "sample_3")

    missing, ploidy = bgen2.read(2,return_probabilities=False,return_missings=True,return_ploidies=True)
    assert_allclose(ploidy[:,0], [1, 2, 2, 2])
    assert_allclose(missing[:,0], [0, 0, 0, 0])
    assert_allclose(bgen2.phased, [0, 1, 1, 0, 1, 1, 1, 1, 0, 0])

if __name__ == '__main__': #!!!cmk remove?
    pytest.main([__file__])