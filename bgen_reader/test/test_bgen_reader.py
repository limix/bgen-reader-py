from __future__ import unicode_literals

import os

import pytest
from bgen_reader import convert_to_dosage, read_bgen, create_metadata_file
from numpy import array
from numpy.testing import assert_, assert_allclose, assert_equal


try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def test_bgen_reader():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    filepath = os.path.join(folder, b"example.32bits.bgen")
    bgen = read_bgen(filepath, verbose=False)
    variants = bgen['variants']
    samples = bgen['samples']
    assert_('genotype' in bgen)

    assert_equal(variants.loc[0, 'chrom'], '01')
    assert_equal(variants.loc[0, 'id'], 'SNPID_2')
    assert_equal(variants.loc[0, 'nalleles'], 2)
    assert_equal(variants.loc[0, 'allele_ids'], "A,G")
    assert_equal(variants.loc[0, 'pos'], 2000)
    assert_equal(variants.loc[0, 'rsid'], 'RSID_2')

    assert_equal(variants.loc[7, 'chrom'], '01')
    assert_equal(variants.loc[7, 'id'], 'SNPID_9')
    assert_equal(variants.loc[7, 'nalleles'], 2)
    assert_equal(variants.loc[7, 'allele_ids'], "A,G")
    assert_equal(variants.loc[7, 'pos'], 9000)
    assert_equal(variants.loc[7, 'rsid'], 'RSID_9')

    n = variants.shape[0]
    assert_equal(variants.loc[n - 1, 'chrom'], '01')
    assert_equal(variants.loc[n - 1, 'id'], 'SNPID_200')
    assert_equal(variants.loc[n - 1, 'nalleles'], 2)
    assert_equal(variants.loc[n - 1, 'allele_ids'], "A,G")
    assert_equal(variants.loc[n - 1, 'pos'], 100001)
    assert_equal(variants.loc[n - 1, 'rsid'], 'RSID_200')

    assert_equal(samples.loc[0, 'id'], 'sample_001')
    assert_equal(samples.loc[7, 'id'], 'sample_008')

    n = samples.shape[0]
    assert_equal(samples.loc[n - 1, 'id'], 'sample_500')

    bgen = read_bgen(filepath, verbose=False)
    variants = bgen['variants']
    assert_('samples' in bgen)
    assert_('genotype' in bgen)

    assert_equal(variants.loc[0, 'chrom'], '01')
    if os.path.exists(filepath + b'.metadata'):
        os.remove(filepath + b'.metadata')


def test_bgen_reader_without_metadata():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    filepath = os.path.join(folder, b"example.32bits.bgen")
    bgen = read_bgen(filepath, verbose=False, metadata_file=False)
    variants = bgen['variants']
    samples = bgen['samples']
    assert_('genotype' in bgen)
    assert_equal(variants.loc[7, 'allele_ids'], "A,G")
    n = samples.shape[0]
    assert_equal(samples.loc[n - 1, 'id'], 'sample_500')
    assert_(not os.path.exists(filepath + b'.metadata'))


def test_bgen_reader_with_wrong_metadata_file():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    filepath = os.path.join(folder, b"example.32bits.bgen")
    metadata_file = os.path.join(folder, b"wrong.metadata")
    with pytest.raises(RuntimeError):
        read_bgen(filepath, verbose=False, metadata_file=metadata_file)


def test_bgen_reader_with_nonexistent_metadata_file():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    filepath = os.path.join(folder, b"example.32bits.bgen")
    metadata_file = os.path.join(folder, b"nonexistent.metadata")
    with pytest.raises(FileNotFoundError):
        read_bgen(filepath, verbose=False, metadata_file=metadata_file)


def test_bgen_reader_file_notfound():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    filepath = os.path.join(folder, b"example.33bits.bgen")
    with pytest.raises(FileNotFoundError):
        read_bgen(filepath, verbose=False)


def test_bgen_reader_convert_to_dosage():
    folder = os.path.dirname(os.path.abspath(__file__)).encode()
    filepath = os.path.join(folder, b"example.32bits.bgen")
    bgen = read_bgen(filepath, verbose=False)
    genotype = bgen['genotype']
    dosage = convert_to_dosage(genotype)
    assert_allclose(dosage[0, 1:3], array([1.93575854, 1.91558579]), rtol=1e-5)
    if os.path.exists(filepath + b'.metadata'):
        os.remove(filepath + b'.metadata')


def test_create_metadata_file():
    folder = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(folder, "example.32bits.bgen")
    metadata_file = os.path.join(folder, filepath + '.metadata')

    create_metadata_file(filepath, metadata_file, verbose=False)
    assert_(os.path.exists(metadata_file))
    os.remove(metadata_file)
