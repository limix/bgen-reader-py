import os
from pathlib import Path

from appdirs import user_cache_dir

from ._file import make_sure_dir_exist

BGEN_READER_CACHE_HOME = Path(
    os.environ.get(
        "BGEN_READER_CACHE_HOME",
        default=Path(user_cache_dir("bgen-reader", "limix")) / "bgen-reader",
    )
)

__all__ = ["BGEN_READER_CACHE_HOME"]

make_sure_dir_exist(BGEN_READER_CACHE_HOME)
make_sure_dir_exist(BGEN_READER_CACHE_HOME / "test_data")
make_sure_dir_exist(BGEN_READER_CACHE_HOME / "metafile")
