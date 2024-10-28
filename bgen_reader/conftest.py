from __future__ import unicode_literals

import pytest


def pytest_configure(config):
    import doctest

    _compatibility()
    pandas_format()
    doctest.ELLIPSIS_MARKER = "-ignore-"

    config.addinivalue_line("markers", "slow: mark test as slow to run")


@pytest.fixture(autouse=True)
def _docdir(request):
    import os

    # Trigger ONLY for the doctests or doctestplus.
    plug = request.config.pluginmanager.getplugin("doctest")
    if plug is None:
        plug = request.config.pluginmanager.getplugin("doctestplus")
        if plug is None:
            item = None
        else:
            item = plug._doctest_textfile_item_cls
    else:
        item = plug.DoctestItem

    if isinstance(request.node, item):
        # Get the fixture dynamically by its name.
        tmpdir = request.getfixturevalue("tmpdir")

        # Chdir only for the duration of the test.
        olddir = os.getcwd()
        tmpdir.chdir()
        yield
        os.chdir(olddir)

    else:
        # For normal tests, we have to yield, since this is a yield-fixture.
        yield


def pandas_format():
    import pandas as pd

    pd.set_option("display.width", 88)
    pd.set_option("display.max_columns", 79)
    pd.set_option("display.max_rows", 8)
    pd.set_option("display.large_repr", "truncate")
    pd.set_option("display.float_format", "{:8.5f}".format)


def _compatibility():
    import warnings

    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


@pytest.fixture
def large_bgen_filepath():
    from pathlib import Path
    from subprocess import check_call

    from bgen_reader._environment import BGEN_READER_CACHE_HOME
    from bgen_reader._file import file_hash

    filepath = BGEN_READER_CACHE_HOME / "test_data" / "large.bgen"

    expected = "b9e75b6c5c5e8c5e1ebd1b2f54a731aef99cfee628895d5008c019716b5909dc"
    if filepath.exists() and file_hash(filepath) == expected:
        return filepath

    pass_filepath = Path("/Users/horta/pass")
    if not pass_filepath.exists():
        raise RuntimeError(f"Could not find {pass_filepath} file.")

    cmd = (
        "curl https://bgen-examples.s3.amazonaws.com/bgen-examples/large.bgen.bz2.enc -s | "
        "openssl enc -d -pbkdf2 -aes-256-cbc -kfile /Users/horta/pass |"
        f"bunzip2 > {filepath}"
    )
    check_call(cmd, shell=True)

    return filepath
