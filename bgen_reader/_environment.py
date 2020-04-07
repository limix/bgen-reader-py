from pathlib import Path

from appdirs import user_cache_dir

from ._file import make_sure_dir_exist

BGEN_READER_CACHE_HOME = Path(user_cache_dir("bgen-reader", "limix")) / "bgen-reader"

breakpoint()
__all__ = ["BGEN_READER_CACHE_HOME"]

make_sure_dir_exist(BGEN_READER_CACHE_HOME)
make_sure_dir_exist(BGEN_READER_CACHE_HOME / "test_data")
make_sure_dir_exist(BGEN_READER_CACHE_HOME / "metafile")
