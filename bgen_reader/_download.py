import os
from urllib.parse import urlparse
from urllib.request import urlretrieve


def download(url):

    dest = os.getcwd()
    _filename = os.path.basename(urlparse(url).path)
    filepath = os.path.join(dest, _filename)
    urlretrieve(url, filepath)
