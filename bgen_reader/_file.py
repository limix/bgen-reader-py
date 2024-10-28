import errno
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path


def make_sure_dir_exist(dirpath: Path):
    dirpath.mkdir(parents=True, exist_ok=True)


def assert_file_exist(filepath: Path):
    if not filepath.exists():
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)


@contextmanager
def tmp_cwd():
    """
    Create and enter a temporary directory.

    The previous working directory is saved and switched back when
    leaving the context. The temporary directory is also recursively
    removed at the context ending.
    """
    oldpwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:

        os.chdir(tmpdir)
        try:
            yield
        finally:
            os.chdir(oldpwd)


def assert_file_readable(filepath: Path):
    with open(filepath, "rb") as f:
        f.read(1)


def is_file_writable(filepath: Path):
    try:
        _touch(filepath)
    except PermissionError:
        return False
    finally:
        if filepath.exists():
            os.remove(filepath)
    return True


def path_to_filename(path: Path):
    drive, rest = path.parts[0], path.parts[1:]
    drive = drive.replace("/", "").replace("\\", "").replace(":", "")
    if len(drive) == 0:
        parts = rest
    else:
        parts = (drive,) + rest
    return Path("%".join(parts))


def file_hash(filepath: Path) -> str:
    import hashlib

    BLOCK_SIZE = 65536

    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        fb = f.read(BLOCK_SIZE)
        while len(fb) > 0:
            sha256.update(fb)
            fb = f.read(BLOCK_SIZE)

    return sha256.hexdigest()


def _touch(filepath: Path, mode=0o666, dir_fd=None, **kwargs):
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
