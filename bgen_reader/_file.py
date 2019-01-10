import errno
import os
import stat
import tempfile
from os.path import exists


def assert_file_exist(filepath):
    if not exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)


def assert_file_readable(filepath):
    with open(filepath, "rb") as f:
        f.read(1)


def permission_write_file(filepath):
    try:
        _touch(filepath)
    except PermissionError:
        return False
    finally:
        if os.path.exists(filepath):
            os.remove(filepath)
    return True


def _get_temp_filepath(folder, filename):
    folder = os.path.abspath(folder)
    folder = os.path.normpath(folder)
    folder_pathname = "%".join(folder.split(os.sep)) + filename
    return os.path.join(tempfile.mkdtemp(), folder_pathname)


def _touch(fname, mode=0o666, dir_fd=None, **kwargs):
    """ Touch a file.

    Credits to <https://stackoverflow.com/a/1160227>.
    """
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(
            f.fileno() if os.utime in os.supports_fd else fname,
            dir_fd=None if os.supports_fd else dir_fd,
            **kwargs,
        )


def _is_group_readable(filepath):
    return bool(os.stat(filepath).st_mode & stat.S_IRGRP)
