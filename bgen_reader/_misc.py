import errno
import os
import stat
from os.path import exists


def _group_readable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IRGRP)


# def make_sure_str(p):
#     try:
#         p = p.decode()
#     except AttributeError:
#         pass
#     return p


def check_file_exist(filepath):
    if not exists(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)


def check_file_readable(filepath):
    if not _group_readable(filepath):
        msg = f"You don't have file permission for reading {filepath}."
        raise RuntimeError(msg)
