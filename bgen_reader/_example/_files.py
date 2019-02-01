import shutil
import tempfile
import warnings
from os.path import dirname, exists, join, realpath
from subprocess import check_call


class example_files(object):
    """ Create a temporary folder with the given files.

    Attributes
    ----------
    filenames : List of available files.
    """

    def __init__(self, filenames):
        self._unlist = False
        if not isinstance(filenames, (tuple, list)):
            filenames = [filenames]
            self._unlist = True

        for fn in filenames:
            if fn not in self.filenames:
                raise ValueError(
                    "Unrecognized file name {}. Choose one of these: {}".format(
                        fn, self.filenames
                    )
                )

        self._dirpath = tempfile.mkdtemp()
        self._filenames = filenames

    filenames = [
        "complex.23bits.bgen",
        "complex.23bits.bgen.metadata",
        "complex.sample",
        "example.32bits.bgen",
        "example.32bits.bgen.metadata",
        "haplotypes.bgen",
        "haplotypes.bgen.metadata.valid",
        "haplotypes.bgen.metadata.corrupted",
        "wrong.metadata",
        "complex.23bits.no.samples.bgen",
        "large.bgen",
    ]

    @property
    def filepath(self):
        return self.__enter__()

    def close(self):
        self.__exit__()

    def __enter__(self):
        import pkg_resources

        filepaths = [join(self._dirpath, fn) for fn in self._filenames]

        for fn, fp in zip(self._filenames, filepaths):
            if fn == "large.bgen":
                cmd = f"cp /Users/horta/code/bgen-speed/large.bgen {fp}"
                cmd = (
                    "curl http://rest.s3for.me/bgen/large.bgen.bz2.enc -s | "
                    "openssl enc -d -pbkdf2 -aes-256-cbc -kfile /Users/horta/pass |"
                    f"bunzip2 > {fp}"
                )
                check_call(cmd, shell=True)
            elif __name__ == "__main__":
                shutil.copy(join(dirname(realpath(__file__)), fn), fp)
            else:
                resource_path = "_example/{}".format(fn)
                content = pkg_resources.resource_string(
                    __name__.split(".")[0], resource_path
                )

                with open(fp, "wb") as f:
                    f.write(content)

        if self._unlist:
            return filepaths[0]
        return filepaths

    def __exit__(self, *_):
        try:
            shutil.rmtree(self._dirpath)
        except PermissionError as e:
            warnings.warn(str(e) + "\n. I will ignore it and proceed.")


def can_run_with(filenames):
    if not isinstance(filenames, (list, tuple)):
        filenames = [filenames]

    if "large.bgen" in filenames:
        ok = shutil.which("curl") is not None and shutil.which("openssl") is not None
        ok = ok and exists("/Users/horta/pass")
        return ok

    return True
