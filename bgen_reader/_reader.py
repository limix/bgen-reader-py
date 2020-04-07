import os
import warnings
from pathlib import Path
from typing import Optional, Union

from ._bgen_file import bgen_file
from ._bgen_metafile import bgen_metafile
from ._environment import BGEN_READER_CACHE_HOME
from ._file import (
    assert_file_exist,
    assert_file_readable,
    is_file_writable,
    path_to_filename,
)
from ._genotype import create_genotypes
from ._metafile import create_metafile
from ._samples import generate_samples, read_samples_file


def read_bgen(
    filepath: Union[str, Path],
    metafile_filepath: Optional[Union[str, Path]] = None,
    samples_filepath: Optional[Union[str, Path]] = None,
    verbose: bool = True,
):
    """
    Read a given BGEN file.

    Parameters
    ----------
    filepath
        Bgen file path.
    metafile_filepath
        File path to the corresponding metafile. A metafile can be created by calling
        :func:`bgen_reader.create_metafile`. If ``None``, a metafile will be automatically created.
        Defaults to ``None``.
    samples_filepath
        Path to a `sample format`_ file or ``None`` to read samples from the bgen file itself.
        Defaults to ``None``.
    verbose
        ``True`` to show progress; ``False`` otherwise. Defaults to ``True``.

    Returns
    -------
    variants : :class:`dask.dataFrame.DataFrame`
        Variant position, chromosomes, rsids, etc.
    samples : :class:`pandas.Series`
        Sample identifications.
    genotype : list
        List of genotypes.

    Examples
    --------
    .. doctest::

        >>> from bgen_reader import example_filepath, read_bgen
        >>>
        >>> bgen = read_bgen(example_filepath("haplotypes.bgen"), verbose=False)
        >>> variants = bgen["variants"]
        >>> samples = bgen["samples"]
        >>>
        >>> v = variants.loc[0].compute()
        >>> g = bgen["genotype"][0].compute()
        >>> print(v)
             id rsid chrom  pos  nalleles allele_ids  vaddr
        0  SNP1  RS1     1    1         2        A,G    102
        >>> print(samples)
        0    sample_0
        1    sample_1
        2    sample_2
        3    sample_3
        Name: id, dtype: object
        >>> print(g["probs"][0])
        [1. 0. 1. 0.]

    .. _sample format: https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html
    """

    filepath = Path(filepath)
    assert_file_exist(filepath)
    assert_file_readable(filepath)

    if metafile_filepath is None:
        metafile_filepath = _infer_metafile_filepath(filepath)
    else:
        metafile_filepath = Path(metafile_filepath)
        assert_file_exist(metafile_filepath)
        assert_file_readable(filepath)

    if not metafile_filepath.exists():
        if verbose:
            print(
                f"We will create the metafile `{metafile_filepath}`. This file will "
                "speed up further\nreads and only need to be created once. So, please, "
                "bear with me."
            )
        create_metafile(filepath, metafile_filepath, verbose)
    elif os.path.getmtime(metafile_filepath) < os.path.getmtime(filepath):
        from ._bgen_metafile import cache as metacache
        from ._genotype import cache as bgencache

        metacache.clear()
        bgencache.clear()

        if verbose:
            print(
                f"File `{filepath}` has been modified after the creation of `{metafile_filepath}`."
                "\nWe will therefore recreate the metadata file. So, please, bear with me."
            )
        os.unlink(metafile_filepath)
        create_metafile(filepath, metafile_filepath, verbose)

    with bgen_file(filepath) as bgen:

        if samples_filepath is None:
            if bgen.contain_samples:
                samples = bgen.read_samples()
            else:
                samples = generate_samples(bgen.nsamples)
        else:
            samples_filepath = Path(samples_filepath)
            assert_file_exist(samples_filepath)
            assert_file_readable(samples_filepath)
            samples = read_samples_file(samples_filepath, verbose)

        with bgen_metafile(metafile_filepath) as metafile:
            variants = metafile.create_variants()

        genotype = create_genotypes(bgen, metafile_filepath, verbose)

    return dict(variants=variants, samples=samples, genotype=genotype)


_metafile_nowrite_dir = """\
You don't have permission to write `{filepath}`. This might prevent speeding-up the reading process
in future runs.
"""


def _infer_metafile_filepath(bgen_filepath: Path) -> Path:
    metafile = bgen_filepath.with_suffix(bgen_filepath.suffix + ".metafile")
    if metafile.exists():
        try:
            assert_file_readable(metafile)
            return metafile
        except RuntimeError as e:
            warnings.warn(str(e), UserWarning)
            return BGEN_READER_CACHE_HOME / "metafile" / path_to_filename(metafile)
    else:
        if is_file_writable(metafile):
            return metafile

        warnings.warn(_metafile_nowrite_dir.format(filepath=metafile), UserWarning)
        return BGEN_READER_CACHE_HOME / "metafile" / path_to_filename(metafile)
