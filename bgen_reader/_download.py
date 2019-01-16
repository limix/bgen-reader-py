import os
import sys
import warnings
from urllib.parse import urlparse
from urllib.request import urlretrieve


def download(url, dest=None, verbose=True, force=False):

    if dest is None:
        dest = os.getcwd()

    _filename = os.path.basename(urlparse(url).path)

    filepath = os.path.join(dest, _filename(url))
    if not force and os.path.exists(filepath):
        warnings.warn("File {} already exists.".format(filepath))
        print("Set `force` to `True` in order to overwrite the existing file.")
        return

    if verbose:
        sys.stdout.write(f"Downloading {url}... ")
    urlretrieve(url, filepath)
    if verbose:
        print("done.")
