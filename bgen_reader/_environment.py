from xdg import XDG_CACHE_HOME

from ._file import make_sure_dir_exist

BGEN_READER_CACHE_HOME = XDG_CACHE_HOME / "bgen-reader"

__all__ = ["BGEN_READER_CACHE_HOME"]

make_sure_dir_exist(BGEN_READER_CACHE_HOME)
make_sure_dir_exist(BGEN_READER_CACHE_HOME / "test_data")
make_sure_dir_exist(BGEN_READER_CACHE_HOME / "metafile")
