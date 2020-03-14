import errno
import os
import stat
import tempfile
from os.path import exists
from pathlib import Path
from xdg import XDG_CACHE_HOME

BGEN_CACHE_HOME = XDG_CACHE_HOME / "bgen"


def make_sure_dir_exist(dirpath: Path):
    dirpath.mkdir(parents=True, exist_ok=True)


make_sure_dir_exist(BGEN_CACHE_HOME)


def assert_file_exist(filepath):
    if not exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)


def assert_file_exist2(filepath: Path):
    if not exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)


def assert_file_readable(filepath):
    with open(filepath, "rb") as f:
        f.read(1)


def assert_file_readable2(filepath: Path):
    with open(filepath, "rb") as f:
        f.read(1)


def is_file_writable(filepath):
    try:
        _touch(filepath)
    except PermissionError:
        return False
    finally:
        if os.path.exists(filepath):
            os.remove(filepath)
    return True


def _get_temp_filepath(path: Path):
    return tempfile.mkdtemp() / Path("%" + "%".join(path.parts))


def path_to_filename(path: Path):
    drive, rest = path.parts[0], path.parts[1:]
    drive = drive.replace("/", "").replace("\\", "").replace(":", "")
    if len(drive) == 0:
        parts = rest
    else:
        parts = (drive,) + rest
    return Path("%".join(parts))


def _touch(filepath, mode=0o666, dir_fd=None, **kwargs):
    """
    Touch a file.

    Credits to <https://stackoverflow.com/a/1160227>.
    """
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(filepath, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(
            f.fileno() if os.utime in os.supports_fd else filepath,
            dir_fd=None if os.supports_fd else dir_fd,
            **kwargs,
        )


def _is_group_readable(filepath):
    return bool(os.stat(filepath).st_mode & stat.S_IRGRP)

