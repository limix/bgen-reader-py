import os
import platform
from shutil import copyfile

import numpy as np
import pytest
from numpy import array, array_equal, isnan
from numpy.testing import assert_allclose, assert_equal

from bgen_reader import example_filepath, open_bgen
from bgen_reader._environment import BGEN_READER_CACHE_HOME
from bgen_reader.test.test_bgen_reader import noread_permission
from bgen_reader.test.write_random import _write_random


def example_filepath2(filename):
    filepath = example_filepath(filename)
    for allow_complex in [False, True]:
        metadata2_path = open_bgen._metadata_path_from_filename(
            filepath, samples_filepath=None, allow_complex=allow_complex
        )
        if metadata2_path.exists():
            metadata2_path.unlink()
    return filepath


def test_bgen_samples_inside_bgen():
    with open_bgen(example_filepath2("haplotypes.bgen"), verbose=False) as _:
        pass


def test_typing():
    with open_bgen(example_filepath2("haplotypes.bgen"), verbose=False) as data:
        with pytest.raises(TypeError):
            data.read(dtype=3)


def test_bgen_samples_not_present():
    with open_bgen(
        example_filepath2("complex.23bits.no.samples.bgen"),
        allow_complex=True,
        verbose=False,
    ) as data:
        samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
        assert all(data.samples == samples)


def test_bgen_samples_specify_samples_file():
    with open_bgen(
        example_filepath2("complex.23bits.bgen"),
        samples_filepath=example_filepath("complex.sample"),
        allow_complex=True,
        verbose=False,
    ) as data:
        samples = ["sample_0", "sample_1", "sample_2", "sample_3"]
        assert all(data.samples == samples)


# TODO: have it back. It was not working anymore.
@pytest.mark.skip
def test_bgen_samples_outside_bgen_unreadable(tmp_path):
    bgen_filepath = example_filepath2("complex.23bits.bgen")
    samples_filepath = tmp_path / "complex.sample"
    copyfile(example_filepath("complex.sample"), samples_filepath)
    with noread_permission(samples_filepath):
        with pytest.raises(PermissionError):
            with open_bgen(
                bgen_filepath, samples_filepath=samples_filepath, verbose=False
            ) as _:
                pass


@pytest.mark.skipif(platform.system() != "Darwin", reason="only reliable on macos")
def test_bgen_file_not_readable(tmp_path):
    filepath = tmp_path / "haplotypes.bgen"
    copyfile(example_filepath2("haplotypes.bgen"), filepath)
    with noread_permission(filepath):
        with pytest.raises(PermissionError):
            with open_bgen(filepath, verbose=False) as _:
                pass


def test_bgen_file_dont_exist():
    with pytest.raises(FileNotFoundError):
        with open_bgen("idontexist.bgen", verbose=False) as _:
            pass


def test_metafile_not_provided():
    with open_bgen(example_filepath2("haplotypes.bgen"), verbose=False) as _:
        pass


def test_open_bgen_phased_genotype():
    filepath = example_filepath2("haplotypes.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:

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

        g = bgen2.read((0, 0))
        assert_allclose(g[0, 0, :], [1.0, 0.0, 1.0, 0.0])
        g = bgen2.read((-1, -1))
        assert_allclose(g[0, 0, :], [1.0, 0.0, 0.0, 1.0])


def test_open_bgen_variants_info():
    filepath = example_filepath2("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:

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

        g = bgen2.read((0, 0))
        assert all(isnan(g[0, 0, :]))

        g = bgen2.read((1, 0))
        a = [0.027802362811705648, 0.00863673794284387, 0.9635608992454505]
        assert_allclose(g[0, 0, :], a)

        b = [
            0.97970582847010945215516,
            0.01947019668749305418287,
            0.00082397484239749366197,
        ]
        g = bgen2.read((2, 1))
        assert_allclose(g[0, 0, :], b)


def test_to_improve_coverage():
    filepath = example_filepath2("example.32bits.bgen")
    bgen2 = open_bgen(filepath, verbose=False)  # Creates metadata2.mmm file
    assert_equal(bgen2.ncombinations[-1], 3)
    assert_equal(bgen2.phased[-1], False)
    with open_bgen(filepath) as bgen2:  # Reuses metadata2.mmm file
        assert_equal(str(bgen2), "open_bgen('{0}')".format(filepath.name))
        assert_equal(bgen2.nsamples, 500)
        assert_equal(bgen2.nvariants, 199)
        assert_equal(bgen2.shape, (500, 199, 3))
        assert_equal(bgen2.ids[-1], "SNPID_200")
        assert_equal(bgen2.rsids[-1], "RSID_200")
        assert_equal(bgen2.chromosomes[-1], "01")
        assert_equal(bgen2.positions[-1], 100001)
        assert_equal(bgen2.nalleles[-1], 2)
        assert_equal(bgen2.allele_ids[-1], "A,G")
        assert_equal(bgen2.ncombinations[-1], 3)
        assert_equal(bgen2.phased[-1], False)
        assert_equal(bgen2.samples[-1], "sample_500")

        b = [
            0.97970582847010945215516,
            0.01947019668749305418287,
            0.00082397484239749366197,
        ]
        g = bgen2.read((2, 1))
        assert_allclose(g[0, 0, :], b)

        g = bgen2.read()
        assert_allclose(g[2, 1, :], b)

    # confirm that out-of-date metadata2 file will be updated
    metadata2 = bgen2._metadata2_path
    del bgen2
    assert os.path.getmtime(metadata2) >= os.path.getmtime(filepath)
    filepath.touch()
    assert os.path.getmtime(metadata2) <= os.path.getmtime(filepath)
    bgen2 = open_bgen(filepath, verbose=False)  # Creates metadata2.mmm file
    del bgen2
    assert os.path.getmtime(metadata2) >= os.path.getmtime(filepath)


def test_to_improve_coverage2():
    filepath = example_filepath2("complex.bgen")
    samplepath = example_filepath2("complex.sample")
    allow_complex = True

    metadata2_path = open_bgen._metadata_path_from_filename(
        filepath, samples_filepath=samplepath, allow_complex=allow_complex
    )

    if metadata2_path.exists():
        metadata2_path.unlink()
    metadata2_temp = metadata2_path.parent / (metadata2_path.name + ".temp")
    metadata2_temp.touch()  # Create an empty .temp file

    bgen2 = open_bgen(
        filepath,
        samples_filepath=samplepath,
        allow_complex=allow_complex,
        verbose=True,
    )  # Creates metadata2.mmm file
    bgen2 = open_bgen(
        filepath,
        samples_filepath=samplepath,
        allow_complex=True,
        verbose=True,
    )  # Creates metadata2.mmm file

    del bgen2


@pytest.mark.skipif(
    "QCTOOLPATH" not in os.environ, reason="This test requires external QCTOOL"
)
@pytest.mark.slow  # It takes hours to generate data locally.  After that, it takes a few minutes
# to run.
def test_bigfile(verbose=False):
    random_file_tests(nsamples=2500, nvariants=500 * 1000, bits=16)


@pytest.mark.skipif(
    "QCTOOLPATH" not in os.environ, reason="This test requires external QCTOOL"
)
@pytest.mark.slow  # Skipping this one by default because it requires the QCTOOL
def test_small_random_file(verbose=False):
    random_file_tests(
        nsamples=25, nvariants=1000, bits=8, verbose=verbose, overwrite=True
    )


def random_file_tests(nsamples, nvariants, bits, verbose=False, overwrite=False):
    test_data_folder = BGEN_READER_CACHE_HOME / "test_data"
    filepath = test_data_folder / "{0}x{1}.{2}bits.bgen".format(
        nsamples, nvariants, bits
    )
    if overwrite or not filepath.exists():
        _write_random(
            filepath,
            nsamples,
            nvariants,
            bits=bits,
            verbose=verbose,
            cleanup_temp_files=True,
        )
    metadata2_path = open_bgen._metadata_path_from_filename(
        filepath, samples_filepath=None, allow_complex=True
    )
    if metadata2_path.exists():
        metadata2_path.unlink()

    with open_bgen(filepath, verbose=verbose) as bgen2:
        assert bgen2.nsamples == nsamples
        assert bgen2.nvariants == nvariants
        val = bgen2.read(-1)
        assert val.shape == (nsamples, 1, 3)
        mean = np.nanmean(val)
        assert mean != mean or (0 <= mean and mean <= 1)


def test_open_bgen_without_metadata():
    filepath = example_filepath2("example.32bits.bgen")
    with open_bgen(filepath, allow_complex=True, verbose=False) as bgen2:
        assert_equal(bgen2.allele_ids[7], "A,G")
        assert_equal(bgen2.samples[-1], "sample_500")


def test_open_bgen_file_notfound():
    with pytest.raises(FileNotFoundError):
        with open_bgen("/1/2/3/example.32bits.bgen", verbose=False) as _:
            pass


def test_open_bgen_complex():
    filepath = example_filepath2("complex.23bits.bgen")
    with open_bgen(filepath, allow_complex=True, verbose=False) as bgen2:

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

        g = bgen2.read((0, 0))
        assert_allclose(g[0, 0, :2], [1, 0])
        assert isnan(g[0, 0, 2])

        g = bgen2.read((1, 0))
        assert_allclose(g[0, 0, :3], [1, 0, 0])

        g = bgen2.read((-1, -1))
        assert_allclose(g[0, 0, :5], [0, 0, 0, 1, 0])

        ploidy = bgen2.read(0, return_probabilities=False, return_ploidies=True)
        assert_allclose(ploidy[:, 0], [1, 2, 2, 2])
        ploidy = bgen2.read(-1, return_probabilities=False, return_ploidies=True)
        assert_allclose(ploidy[:, 0], [4, 4, 4, 4])

        assert_equal(bgen2.phased.dtype.name, "bool")
        ideal = array([False, True, True, False, True, True, True, True, False, False])
        assert array_equal(bgen2.phased, ideal)


def test_open_bgen_complex_sample_file():
    with open_bgen(
        example_filepath2("complex.23bits.bgen"),
        samples_filepath=example_filepath("complex.sample"),
        allow_complex=True,
        verbose=False,
    ) as bgen2:

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

        missing, ploidy = bgen2.read(
            2, return_probabilities=False, return_missings=True, return_ploidies=True
        )
        assert_allclose(ploidy[:, 0], [1, 2, 2, 2])
        assert_allclose(missing[:, 0], [0, 0, 0, 0])
        assert_allclose(bgen2.phased, [0, 1, 1, 0, 1, 1, 1, 1, 0, 0])


def test_close_del_with():
    filepath = example_filepath2("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        pass
    with pytest.raises(ValueError):
        bgen2.read()

    bgen2 = open_bgen(filepath, verbose=False)
    bgen2.close()
    with pytest.raises(ValueError):
        bgen2.read()


def test_read_max_combinations():
    filepath = example_filepath2("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        assert (
            np.mean(np.isnan(bgen2.read())) < 2.1e-05
        )  # main data as only a few missing
        val = bgen2.read(max_combinations=5)
        assert (
            np.mean(np.isnan(val[:, :, :3])) < 2.1e-05
        )  # main data as only a few missing
        assert np.all(np.isnan(val[:, :, 3:]))  # all the extra are NaN
        with pytest.raises(ValueError):
            bgen2.read(max_combinations=2)
        with pytest.raises(ValueError):
            bgen2.read(max_combinations=1)
        with pytest.raises(ValueError):
            bgen2.read(max_combinations=0)


def test_read_dtype_and_order():
    filepath = example_filepath2("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        full = bgen2.read()
        assert full.dtype == np.float64
        assert full.flags["F_CONTIGUOUS"] and not full.flags["C_CONTIGUOUS"]

        val = bgen2.read(None, dtype="float32", order="C")
        assert val.dtype == np.float32
        assert val.flags["C_CONTIGUOUS"] and not val.flags["F_CONTIGUOUS"]
        assert np.allclose(full, val, atol=5e-8, equal_nan=True)


def test_read_indexing():
    filepath = example_filepath2("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        full = bgen2.read()

        val = bgen2.read(22)
        assert np.allclose(full[:, [22]], val, equal_nan=True)

        val = bgen2.read([22])
        assert np.allclose(full[:, [22]], val, equal_nan=True)

        val = bgen2.read([22, 30])
        assert np.allclose(full[:, [22, 30]], val, equal_nan=True)

        val = bgen2.read(slice(10, 30, 2))
        assert np.allclose(full[:, 10:30:2], val, equal_nan=True)

        bool_list = [i % 2 == 0 for i in range(bgen2.nvariants)]
        val = bgen2.read(bool_list)
        assert np.allclose(full[:, bool_list], val, equal_nan=True)

        val = bgen2.read((None, None))
        assert np.allclose(full, val, equal_nan=True)

        val = bgen2.read((22, None))
        assert np.allclose(full[[22], :], val, equal_nan=True)

        val = bgen2.read((22, [11, 9]))
        assert np.allclose(full[[22], :][:, [11, 9]], val, equal_nan=True)

        val = bgen2.read(([22, 30], [11, 9]))
        assert np.allclose(full[[22, 30], :][:, [11, 9]], val, equal_nan=True)

        val = bgen2.read((slice(10, 30, 2), [11, 9]))
        assert np.allclose(full[10:30:2, :][:, [11, 9]], val, equal_nan=True)

        bool_list = [i % 2 == 0 for i in range(bgen2.nsamples)]
        val = bgen2.read((bool_list, [11, 9]))
        assert np.allclose(full[bool_list, :][:, [11, 9]], val, equal_nan=True)

        val = bgen2.read(([-1], [-1]))
        assert np.allclose(full[-1, -1], val, equal_nan=True)

        val = bgen2.read(np.s_[10:30:2, :5])
        assert np.allclose(full[10:30:2, :5, :], val, equal_nan=True)

        # Read no variants
        val, missing, ploidy = bgen2.read(
            [], return_missings=True, return_ploidies=True
        )
        assert val.shape == (bgen2.nsamples, 0, bgen2.max_combinations)
        assert missing.shape == (bgen2.nsamples, 0)
        assert ploidy.shape == (bgen2.nsamples, 0)

        # Read no samples and no variants
        val, missing, ploidy = bgen2.read(
            ([], []), return_missings=True, return_ploidies=True
        )
        assert val.shape == (0, 0, bgen2.max_combinations)
        assert missing.shape == (0, 0)
        assert ploidy.shape == (0, 0)


def test_read_multiple_returns():
    filepath = example_filepath2("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        full, full_missing, full_ploidy = bgen2.read(
            return_missings=True, return_ploidies=True
        )

        val, missing = bgen2.read(return_missings=True)
        assert np.allclose(full, val, equal_nan=True)
        assert np.allclose(full_missing, missing, equal_nan=False)

        ploidy = bgen2.read(return_probabilities=False, return_ploidies=True)
        assert np.allclose(full_ploidy, ploidy, equal_nan=False)

        val, missing = bgen2.read((slice(10, 30, 2), [11, 9]), return_missings=True)
        assert np.allclose(full[10:30:2, :][:, [11, 9]], val, equal_nan=True)
        assert np.allclose(
            full_missing[10:30:2, :][:, [11, 9]], missing, equal_nan=False
        )

        ploidy = bgen2.read(
            (slice(10, 30, 2), [11, 9]),
            return_probabilities=False,
            return_ploidies=True,
        )
        assert np.allclose(full_ploidy[10:30:2, :][:, [11, 9]], ploidy, equal_nan=False)


def test_coverage3():
    with pytest.raises(ValueError):
        with open_bgen(
            example_filepath2("example.bgen"),
            samples_filepath=example_filepath(
                "complex.sample"
            ),  # Wrong size sample file
            verbose=False,
        ) as _:
            pass

    with pytest.raises(ValueError):
        with open_bgen(
            example_filepath2("complex.bgen"),
            verbose=False,
        ) as _:
            pass


def test_coverage4(tmp_path):
    oldpwd = os.getcwd()
    filepath = example_filepath2("example.32bits.bgen")
    try:
        os.chdir(filepath.parent)
        with open_bgen(filepath.name) as bgen2:
            assert bgen2.shape == (500, 199, 3)
    finally:
        os.chdir(oldpwd)


def test_allele_expectation():
    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        e = bgen2.allele_expectation(
            np.s_[bgen2.samples == "sample_005", bgen2.rsids == "RSID_6"]
        )
        assert np.allclose(e, [[[1.01086423, 0.98913577]]])

    with pytest.raises(ValueError):
        filepath = example_filepath("haplotypes.bgen")
        with open_bgen(filepath, verbose=False) as bgen2:
            bgen2.allele_expectation()

    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        e = bgen2.allele_expectation(np.s_[:, []])
        assert e.shape == (500, 0, 2)

    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        e = bgen2.allele_expectation(
            np.s_[bgen2.samples == "sample_005", bgen2.rsids == "RSID_6"],
            assume_constant_ploidy=False,
        )
        assert np.allclose(e, [[[1.01086423, 0.98913577]]])

    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        e = bgen2.allele_expectation(np.s_[:, []], assume_constant_ploidy=False)
        assert e.shape == (500, 0, 2)


def test_threads():
    filepath = example_filepath("example.32bits.bgen")
    with open_bgen(filepath, verbose=False) as bgen2:
        for num_threads in [1, 2]:
            for slice in [np.s_[:, :], np.s_[:, []]]:
                val = bgen2.read(index=slice, num_threads=num_threads)
                row_count = len(bgen2.samples[slice[0]])
                col_count = len(bgen2.ids[slice[1]])
                assert val.shape == (row_count, col_count, 3)


# def test_move_mmm():
#     filepath = example_filepath("example.32bits.bgen")
#     with open_bgen(
#         filepath, metadata_filepath=r"m:\deldir\cmk.mmm", verbose=False
#     ) as bgen2:
#         for num_threads in [1, 2]:
#             for slice in [np.s_[:, :], np.s_[:, []]]:
#                 val = bgen2.read(index=slice, num_threads=num_threads)
#                 row_count = len(bgen2.samples[slice[0]])
#                 col_count = len(bgen2.ids[slice[1]])
#                 assert val.shape == (row_count, col_count, 3)
#     with open_bgen(
#         filepath, metadata_filepath=r"m:\deldir\cmk.mmm", verbose=False
#     ) as bgen2:
#         for num_threads in [1, 2]:
#             for slice in [np.s_[:, :], np.s_[:, []]]:
#                 val = bgen2.read(index=slice, num_threads=num_threads)
#                 row_count = len(bgen2.samples[slice[0]])
#                 col_count = len(bgen2.ids[slice[1]])
#                 assert val.shape == (row_count, col_count, 3)


if __name__ == "__main__":
    pytest.main([__file__])
