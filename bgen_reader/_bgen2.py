import hashlib
import math
import multiprocessing
import os
import threading
from os.path import getmtime
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import numpy as np
from cbgen import bgen_file, bgen_metafile

# from ._bgen_metafile import bgen_metafile
from cbgen._ffi import ffi as cb_ffi
from cbgen._ffi import lib as cb_lib
from numpy import asarray, stack

from ._file import assert_file_exist, assert_file_readable, tmp_cwd
from ._helper import _log_in_place, genotypes_to_allele_counts, get_genotypes
from ._metafile import infer_metafile_filepath
from ._multimemmap import MultiMemMap


# https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
class open_bgen:
    """
    A NumPy-inspired class for fast opening and reading of BGEN files.

    Parameters
    ----------
    filepath
        BGEN file path.
    samples_filepath
        Path to a `sample format`_ file or ``None`` to read samples from the BGEN file itself.
        Defaults to ``None``.
    allow_complex
        ``False`` (default) to assume homogeneous data; ``True`` to allow complex data.
        The BGEN format allows every variant to vary in its phaseness, its allele count,
        and its maximum ploidy. For files where these values may actually vary,
        set ``allow_complex`` to ``True``.

    verbose
        ``True`` (default) to show progress; ``False`` otherwise.

    Returns
    -------
    an open_bgen object : :class:`open_bgen`


    The first time a file is opened , ``open_bgen`` creates a .metadata2.mmm file, a process that takes seconds to hours,
    depending on the size of the file and the ``allow_complex`` setting. Subsequent openings take just a fraction of
    a second. Changing ``samples_filepath`` or ``allow_complex`` results in a new .metadata2.mmm with a slightly
    different name.

    .. _open_examples:

    Examples
    --------
    With the `with <https://docs.python.org/3/reference/compound_stmts.html#grammar-token-with-stmt>`__ statement, list :attr:`samples` and variant :attr:`ids`, then :meth:`read` the whole file.

    .. doctest::

        >>> from bgen_reader import example_filepath, open_bgen
        >>>
        >>> file = example_filepath("haplotypes.bgen")
        >>> with open_bgen(file, verbose=False) as bgen:
        ...     print(bgen.ids)
        ...     print(bgen.samples)
        ...     print(bgen.read())
        ['SNP1' 'SNP2' 'SNP3' 'SNP4']
        ['sample_0' 'sample_1' 'sample_2' 'sample_3']
        [[[1. 0. 1. 0.]
          [0. 1. 1. 0.]
          [1. 0. 0. 1.]
          [0. 1. 0. 1.]]
        <BLANKLINE>
         [[0. 1. 1. 0.]
          [1. 0. 0. 1.]
          [0. 1. 0. 1.]
          [1. 0. 1. 0.]]
        <BLANKLINE>
         [[1. 0. 0. 1.]
          [0. 1. 0. 1.]
          [1. 0. 1. 0.]
          [0. 1. 1. 0.]]
        <BLANKLINE>
         [[0. 1. 0. 1.]
          [1. 0. 1. 0.]
          [0. 1. 1. 0.]
          [1. 0. 0. 1.]]]

    Open the file (without `with`) and read probabilities for one variant.

    .. doctest::

        >>> bgen = open_bgen(file, verbose=False)
        >>> print(bgen.read(2))
        [[[1. 0. 0. 1.]]
        <BLANKLINE>
         [[0. 1. 0. 1.]]
        <BLANKLINE>
         [[1. 0. 1. 0.]]
        <BLANKLINE>
         [[0. 1. 1. 0.]]]
        >>> del bgen                 # close and delete object

    Open the file and then first read for a :class:`slice` of samples and variants, and then for a single sample and variant.

    .. doctest::

        >>> bgen = open_bgen(file, verbose=False)
        >>> print(bgen.read((slice(1,3),slice(2,4))))
        [[[0. 1. 0. 1.]
          [1. 0. 1. 0.]]
        <BLANKLINE>
         [[1. 0. 1. 0.]
          [0. 1. 1. 0.]]]
        >>> print(bgen.read((0,1)))
        [[[0. 1. 1. 0.]]]
        >>> del bgen                 # close and delete object

    .. _sample format: https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        samples_filepath: Optional[Union[str, Path]] = None,
        allow_complex: bool = False,
        verbose: bool = True,
    ):
        filepath = Path(filepath)
        assert_file_exist(filepath)
        assert_file_readable(filepath)
        self._filepath = filepath
        self._allow_complex = allow_complex
        self._verbose = verbose
        self._samples_filepath = (
            Path(samples_filepath) if samples_filepath is not None else None
        )

        # Getting absolute paths is hard: https://discuss.python.org/t/pathlib-absolute-vs-resolve/2573/10
        self._metadata2_path = self._metadata_path_from_filename(
            self._filepath, self._samples_filepath, self._allow_complex
        ).resolve()  # needed because of tmp_cwd in create_metadata
        if not self._metadata2_path.is_absolute():
            self._metadata2_path = Path.cwd() / self._metadata2_path

        self._cbgen = bgen_file(filepath)

        self._nvariants = self._cbgen.nvariants
        self._nsamples = self._cbgen.nsamples

        if self._metadata2_path.exists() and getmtime(self._metadata2_path) < getmtime(
            self._filepath
        ):
            self._metadata2_path.unlink()
        if not self._metadata2_path.exists():
            self._create_metadata2()

        self._metadata2_memmaps = MultiMemMap(self._metadata2_path, mode="r")

    def _create_metadata2(self):
        metadata2_temp = self._metadata2_path.parent / (
            self._metadata2_path.name + ".temp"
        )
        if metadata2_temp.exists():
            metadata2_temp.unlink()

        with MultiMemMap(metadata2_temp, mode="w+") as mmm_wplus:
            self._extract_nalleles_ids_etc(mmm_wplus)
            self._extract_ncombinations_etc(mmm_wplus)
            self._extract_samples_etc(mmm_wplus)
        os.rename(metadata2_temp, self._metadata2_path)

    def _extract_samples_etc(self, mmm_wplus):
        if self._samples_filepath is not None:
            self._extract_samples_from_samples_file(mmm_wplus)
        else:
            if self._cbgen.contain_samples:
                self._extract_samples_from_bgen_file(mmm_wplus)
            else:
                self._extract_samples_from_nothing(mmm_wplus)
        self._extract_sample_range(mmm_wplus)

    def _extract_samples_from_nothing(self, mmm_wplus):
        with _log_in_place("metadata", self._verbose) as updater:
            prefix = "sample_"
            max_length = len(prefix + str(self.nsamples - 1))
            samples_memmap = mmm_wplus.append_empty(
                "samples", self.nsamples, f"<U{max_length}"
            )
            for i in range(
                self.nsamples
            ):  # LATER Is there another low-memory way to do this that would be faster?
                if i % 1000 == 0 or i + 1 == self.nsamples:
                    updater(
                        "'generate samples': part {0:,} of {1:,}".format(
                            i + 1, self.nsamples
                        )
                    )
                samples_memmap[i] = prefix + str(i)

    def _extract_samples_from_bgen_file(self, mmm_wplus):
        with _log_in_place("metadata", self._verbose) as updater:
            # Use cb_lib instead of Python cbgen so can allocate to memory map instead of RAM.
            bgen_samples = cb_lib.bgen_file_read_samples(self._cbgen._bgen_file)
            if bgen_samples == cb_ffi.NULL:
                raise RuntimeError("Could not fetch samples from the bgen file.")

            try:
                samples_max_len = cb_ffi.new("uint32_t[]", 1)
                cb_lib.read_samples_part1(bgen_samples, self.nsamples, samples_max_len)
                updater("'samples from bgen'")
                samples_memmap = mmm_wplus.append_empty(
                    "samples", self.nsamples, f"<U{samples_max_len[0]}"
                )
                samples = mmm_wplus.append_empty(
                    "_samples", self.nsamples, f"S{samples_max_len[0]}"
                )  # This one second, so we can delete it afterwards.
                cb_lib.read_samples_part2(
                    bgen_samples,
                    self.nsamples,
                    cb_ffi.from_buffer("char[]", samples),
                    samples_max_len[0],
                )
            finally:
                cb_lib.bgen_samples_destroy(bgen_samples)
            samples_memmap[:] = samples
            mmm_wplus.popitem()  # Remove _samples

    def _extract_sample_range(self, mmm_wplus):
        sample_range_memmap = mmm_wplus.append_empty(
            "sample_range", self.nsamples, "int32"
        )
        with _log_in_place("metadata", self._verbose) as updater:
            for i in range(
                self.nsamples
            ):  # LATER: Is there another low memory way to do this that would be faster?
                if i % 1000 == 0 or i + 1 == self.nsamples:
                    updater(
                        "'sample range': part {0:,} of {1:,}".format(
                            i + 1, self.nsamples
                        )
                    )
                sample_range_memmap[i] = i

    def _extract_samples_from_samples_file(self, mmm_wplus):
        with _log_in_place("metadata", self._verbose) as updater:
            assert_file_exist(self._samples_filepath)
            assert_file_readable(self._samples_filepath)
            if self._verbose:
                print(f"Sample IDs are read from '{self._samples_filepath}''.")

            # Find max length
            max_len = 0
            with self._samples_filepath.open("r") as fp:
                fp.readline()
                fp.readline()
                for index, line in enumerate(fp):
                    if index % 1000 == 0 or index + 1 == self.nsamples:
                        updater(
                            "'samples_filepath max_len': part {0:,} of {1:,}".format(
                                index + 1, self.nsamples
                            )
                        )
                    max_len = max(max_len, len(line.strip()))

                if index + 1 != self.nsamples:
                    raise ValueError(
                        "Expect number of samples in file to match number of samples in BGEN file"
                    )

            # Copy samples into memmap
            samples_memmap = mmm_wplus.append_empty(
                "samples", self.nsamples, f"<U{max_len}"
            )
            with self._samples_filepath.open("r") as fp:
                fp.readline()
                fp.readline()
                for index, line in enumerate(fp):
                    if index % 1000 == 0 or index + 1 == self.nsamples:
                        updater(
                            "'samples_filepath': part {0:,} of {1:,}".format(
                                index + 1, self.nsamples
                            )
                        )
                    samples_memmap[index] = line.strip()

    def _extract_ncombinations_etc(self, mmm_wplus):

        ncombinations_memmap = mmm_wplus.append_empty(
            "ncombinations", (self.nvariants), "int32"
        )
        phased_memmap = mmm_wplus.append_empty("phased", (self.nvariants), "bool")

        if not self._allow_complex:
            if self.nvariants == 0:
                raise ValueError(
                    "If allow_complex is False, there must be at least one variant"
                )
            i = 0
            previous_i = None
            while i < self.nvariants:
                genotype = self._cbgen.read_genotype(mmm_wplus["vaddr"][i])
                ncombinations_memmap[i] = genotype.probability.shape[-1]
                phased_memmap[i] = genotype.phased
                if previous_i is not None and (
                    ncombinations_memmap[previous_i] != ncombinations_memmap[i]
                    or phased_memmap[previous_i] != phased_memmap[i]
                ):
                    raise ValueError(
                        "allow_complex is False but a spot check shows that file is complex"
                    )
                previous_i = i
                i = (i + 1) * 2 - 1
            ncombinations_memmap[1:] = ncombinations_memmap[0]
            phased_memmap[1:] = phased_memmap[0]
        else:
            if self._verbose:
                print(
                    "Parameter 'allow_complex' is True, so reading phase and ncombinations of every variant"
                )
            with _log_in_place("metadata", self._verbose) as updater:
                for i, vaddr0 in enumerate(mmm_wplus["vaddr"]):
                    if i % 1000 == 0 or i + 1 == self.nsamples:
                        updater(
                            "'ncombinations': part {0:,} of {1:,}".format(
                                i + 1, self.nvariants
                            )
                        )
                    genotype = None
                    try:
                        genotype = cb_lib.bgen_file_open_genotype(
                            self._cbgen._bgen_file, vaddr0
                        )
                        ncombinations_memmap[i] = cb_lib.bgen_genotype_ncombs(genotype)
                        phased_memmap[i] = cb_lib.bgen_genotype_phased(genotype)
                    finally:
                        if genotype is not None:
                            cb_lib.bgen_genotype_close(genotype)

        max_combinations_memmap = mmm_wplus.append_empty(
            "max_combinations", (1), "int32"
        )
        max_combinations_memmap[0] = max(mmm_wplus["ncombinations"])

    def read(
        self,
        index: Optional[Any] = None,
        dtype: Optional[Union[type, str]] = np.float64,
        order: Optional[str] = "F",
        max_combinations: Optional[int] = None,
        return_probabilities: Optional[bool] = True,
        return_missings: Optional[bool] = False,
        return_ploidies: Optional[bool] = False,
        num_threads: Optional[int] = None,
    ) -> Union[
        None,
        np.ndarray,
        Tuple[np.ndarray, np.ndarray],
        Tuple[np.ndarray, np.ndarray, np.ndarray],
    ]:
        """
        Read genotype information from an :class:`open_bgen` object.

        Parameters
        ----------
        index
            An expression specifying the samples and variants to read. (See :ref:`read_examples`, below).
            Defaults to ``None``, meaning read all.
        dtype : data-type
            The desired data-type for the returned probability array.
            Defaults to :class:`numpy.float64`. Use :class:`numpy.float32` or :class:`numpy.float16`, when appropriate,
            to save 50% or 75% of memory. (See :ref:`read_notes`, below).
        order : {'F','C'}
            The desired memory layout for the returned probability array.
            Defaults to ``F`` (Fortran order, which is variant-major).
        max_combinations : int or ``None``.
            The number of values to allocate for each probability distribution.
            Defaults to a number just large enough for any data in the file.
            For unphased, diploid, biallelic data, it will default to 3. For phased, diploid, biallelic data, it will
            default to 4. Any overallocated space is filled with :const:`numpy.nan`.
        return_probabilities: bool
            Read and return the probabilities for samples and variants specified.
            Defaults to ``True``.
        return_missings: bool
            Return a boolean array telling which probabilities are missing.
            Defaults to ``False``.
        return_ploidies: bool
            Read and return the ploidy for the samples and variants specified.
            Defaults to ``False``.
        num_threads: bool
            The number of threads with which to read data. Defaults to all available processors or the number of variants
            being read, whichever is less. Can also be set with the 'MKL_NUM_THREADS' environment variable.

        Returns
        -------
        zero to three :class:`numpy.ndarray`
            always in this order:

            * a :class:`numpy.ndarray` of probabilities with ``dtype`` and shape `(nsamples_out,nvariants_out,max_combinations)`,
              if ``return_probabilities`` is ``True`` (the default). Missing data is filled with :const:`numpy.nan`.
            * a :class:`numpy.ndarray` of ``bool`` of shape `(nsamples_out,nvariants_out)`, if ``return_missings`` is ``True``
            * a :class:`numpy.ndarray` of ``int`` of shape `(nsamples_out,nvariants_out)`, if ``return_ploidies`` is ``True``


        .. _read_notes:

        Notes
        ------
        * About ``dtype``

            If you know the compression level of your BGEN file, you can sometimes save 50% or 75% on memory with ``dtype``.
            (Test with your data to confirm you are not losing any precision.) The approximate relationship is:

                * BGEN compression 1 to 10 bits: ``dtype`` ='float16'
                * BGEN compression 11 to 23 bits: ``dtype`` ='float32'
                * BGEN compression 24 to 32 bits: ``dtype`` ='float64' (default)


        .. _read_examples:

        Examples
        --------
        * Index Examples

            To read all data in a BGEN file, set ``index`` to ``None``. This is the default.

            .. doctest::

                >>> import numpy as np
                >>> from bgen_reader import example_filepath, open_bgen
                >>>
                >>> with open_bgen(example_filepath("haplotypes.bgen"), verbose=False) as bgen_h:
                ...     print(bgen_h.read()) #real all
                [[[1. 0. 1. 0.]
                  [0. 1. 1. 0.]
                  [1. 0. 0. 1.]
                  [0. 1. 0. 1.]]
                <BLANKLINE>
                 [[0. 1. 1. 0.]
                  [1. 0. 0. 1.]
                  [0. 1. 0. 1.]
                  [1. 0. 1. 0.]]
                <BLANKLINE>
                 [[1. 0. 0. 1.]
                  [0. 1. 0. 1.]
                  [1. 0. 1. 0.]
                  [0. 1. 1. 0.]]
                <BLANKLINE>
                 [[0. 1. 0. 1.]
                  [1. 0. 1. 0.]
                  [0. 1. 1. 0.]
                  [1. 0. 0. 1.]]]

            To read selected variants, set ``index`` to an ``int``, a list of ``int``, a :class:`slice`, or a list of ``bool``.
            Negative integers count from the end of the data.

            .. doctest::

                >>> bgen_e = open_bgen(example_filepath("example.bgen"), verbose=False)
                >>> probs = bgen_e.read(5)  # read the variant indexed by 5.
                >>> print(probs.shape)      # print the dimensions of the returned numpy array.
                (500, 1, 3)
                >>> probs = bgen_e.read([5,6,1])  # read the variant indexed by 5, 6, and 1
                >>> print(probs.shape)
                (500, 3, 3)
                >>> probs = bgen_e.read(slice(5)) #read the first 5 variants
                >>> print(probs.shape)
                (500, 5, 3)
                >>> probs = bgen_e.read(slice(2,5)) #read variants from 2 (inclusive) to 5 (exclusive)
                >>> print(probs.shape)
                (500, 3, 3)
                >>> probs = bgen_e.read(slice(2,None)) # read variants starting at index 2.
                >>> print(probs.shape)
                (500, 197, 3)
                >>> probs = bgen_e.read(slice(None,None,10)) #read every 10th variant
                >>> print(probs.shape)
                (500, 20, 3)
                >>> print(np.unique(bgen_e.chromosomes)) # print unique chrom values
                ['01']
                >>> probs = bgen_e.read(bgen_e.chromosomes=='01') # read all variants in chrom 1
                >>> print(probs.shape)
                (500, 199, 3)
                >>> probs = bgen_e.read(-1) # read the last variant
                >>> print(probs.shape)
                (500, 1, 3)

            To read selected samples, set ``index`` to a tuple of the form ``(sample_index,None)``, where ``sample index`` follows the form
            of ``variant index``, above.

            .. doctest::

                >>> probs = bgen_e.read((0,None)) # Read 1st sample (across all variants)
                >>> print(probs.shape)
                (1, 199, 3)
                >>> probs = bgen_e.read((slice(None,None,10),None)) # Read every 10th sample
                >>> print(probs.shape)
                (50, 199, 3)

            To read selected samples and selected variants, set ``index`` to a tuple of the form ``(sample_index,variant_index)``,
            where ``sample_index`` and ``variant_index`` follow the forms above.

            .. doctest::

                >>> # Read samples 10 (inclusive) to 20 (exclusive) and the first 15 variants.
                >>> probs = bgen_e.read((slice(10,20),slice(15)))
                >>> print(probs.shape)
                (10, 15, 3)
                >>> #read last and 2nd-to-last sample and the last variant
                >>> probs = bgen_e.read(([-1,-2],-1))
                >>> print(probs.shape)
                (2, 1, 3)

        * Multiple Return Example

            Read probabilities, missingness, and ploidy. Print all unique ploidies values.

            .. doctest::

                >>> probs,missing,ploidy = bgen_e.read(return_missings=True,return_ploidies=True)
                >>> print(np.unique(ploidy))
                [2]
        """
        # LATER could allow strings (variant names) and lists of strings
        if not hasattr(self, "_cbgen"):
            raise ValueError("I/O operation on a closed file")

        max_combinations = (
            max_combinations if max_combinations is not None else self.max_combinations
        )  # Can't use 'or' because it treats 0 as False

        samples_index, variants_index = self._split_index(index)

        samples_index = self._sample_range[
            samples_index
        ]  # converts slice(), etc to a list of  numbers
        vaddr = self._vaddr[variants_index]
        num_threads = self._get_num_threads(num_threads, len(vaddr))
        ncombinations = self._ncombinations[variants_index]

        if len(ncombinations) > 0 and max(ncombinations) > max_combinations:
            raise ValueError(
                "Need at least {0} max_combinations, but only {1} given".format(
                    max(ncombinations), max_combinations
                )
            )

        dtype = np.dtype(dtype)
        if return_probabilities:
            val = np.full(
                (len(samples_index), len(vaddr), max_combinations),
                np.nan,
                dtype=dtype,
                order=order,
            )
        if return_missings:
            missing_val = np.full(
                (len(samples_index), len(vaddr)), False, dtype="bool", order=order
            )
        if return_ploidies:
            ploidy_val = np.full(
                (len(samples_index), len(vaddr)), 0, dtype="int", order=order
            )

        approx_read_seconds = len(vaddr) / 20000.0 + len(vaddr) * self.nsamples / (
            2 * 1000 * 1000.0
        )
        vaddr_per_second = max(1, len(vaddr) // int(max(1, approx_read_seconds)))
        vaddr_per_second = 10 ** (
            int(math.log10(vaddr_per_second) + 0.5)
        )  # Do "logarithmic rounding" to make numbers look nicer, e.g.  999 -> 1000
        with _log_in_place("reading", self._verbose) as updater:

            def worker(start, end, thread_index, thread_count):
                with bgen_file(self._filepath) as cbgen:
                    for out_index in range(start, end):
                        vaddr0 = vaddr[out_index]
                        if (
                            thread_index == 0
                            and (out_index + 1) % vaddr_per_second == 0
                            or out_index + 1 == end - start
                        ):
                            updater(
                                f"thread {thread_index+1} of {thread_count}, part {out_index + 1:,} of {end-start:,}"
                            )

                        if dtype == np.float16 or dtype == np.float32:
                            precision = 32
                        else:
                            precision = 64

                        if return_missings or return_ploidies:
                            genotype = cbgen.read_genotype(vaddr0, precision)
                            probability = genotype.probability
                        elif return_probabilities:
                            probability = cbgen.read_probability(vaddr0, precision)

                        if return_probabilities:
                            val[:, out_index, : ncombinations[out_index]] = (
                                probability
                                if (samples_index is self._sample_range)
                                else probability[samples_index, :]
                            )

                        if return_missings:
                            missing_val[:, out_index] = (
                                genotype.missings
                                if (samples_index is self._sample_range)
                                else genotype.missing[samples_index]
                            )

                        if return_ploidies:
                            ploidy_val[:, out_index] = (
                                genotype.ploidy
                                if (samples_index is self._sample_range)
                                else genotype.ploidy[samples_index]
                            )

            threads = []
            vaddr_per_thread = -(-len(vaddr) // num_threads)  # Int Ceiling
            start = 0
            for thread_index in range(num_threads):
                end = min(start + vaddr_per_thread, len(vaddr))
                threads.append(
                    threading.Thread(
                        target=worker, args=(start, end, thread_index, num_threads)
                    )
                )
                start = end
            for t in threads:
                t.start()
            # Wait for all threads to complete
            for t in threads:
                t.join()

        result_array = (
            ([val] if return_probabilities else [])
            + ([missing_val] if return_missings else [])
            + ([ploidy_val] if return_ploidies else [])
        )
        if len(result_array) == 1:
            return result_array[0]
        else:
            return tuple(result_array)

    def _get_num_threads(self, num_threads, len_vaddr):
        if num_threads is not None:
            return max(min(num_threads, len_vaddr), 1)
        if "MKL_NUM_THREADS" in os.environ:
            return max(min(int(os.environ["MKL_NUM_THREADS"]), len_vaddr), 1)
        return max(min(multiprocessing.cpu_count(), len_vaddr), 1)

    @property
    def nsamples(self) -> int:
        """
        The number of samples in the data (``int``).

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.nsamples)
            4

        """
        return self._nsamples

    @property
    def nvariants(self) -> int:
        """
        The number of variants in the data (``int``).

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.nvariants)
            4

        """
        return self._nvariants

    @property
    def max_combinations(self) -> int:
        """
        The maximum number of values in any variant's probability distribution (``int``).

        For unphased, diploidy, biallelic data, it will be 3. For phased, diploidy, biallelic data it will be 4. In general,
        it is the maximum value in :attr:`~bgen_reader.open_bgen.ncombinations`.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.max_combinations)
            4

        """
        return self._metadata2_memmaps["max_combinations"][0]

    @property
    def shape(self) -> Tuple[int, int, int]:
        """
        The tuple (:attr:`~bgen_reader.open_bgen.nsamples`, :attr:`~bgen_reader.open_bgen.nvariants`,
        :attr:`~bgen_reader.open_bgen.max_combinations`).

        Example
        --------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.shape)
            (4, 4, 4)

        """
        return (self.nsamples, self.nvariants, self.max_combinations)

    # This is static so that test code can use it easily.
    # LATER could make a version of this method public
    @staticmethod
    def _metadata_path_from_filename(filename, samples_filepath, allow_complex=False):
        # If there is a sample file, put a hash its name is the name of the metadata file
        if samples_filepath is None:
            s_string = ""
        else:
            hash = hashlib.sha256(samples_filepath.name.encode("utf-8")).hexdigest()
            s_string = ".S" + hash[:6]
        a_string = "" if not allow_complex else ".complex"
        return infer_metafile_filepath(
            Path(filename), f"{s_string}{a_string}.metadata2.mmm"
        )

    @property
    def samples(self) -> List[str]:
        """
        The sample identifiers (a :class:`numpy.ndarray` of ``str``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.samples)
            ['sample_0' 'sample_1' 'sample_2' 'sample_3']

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     samples = bgen.samples.copy()
            >>> print(samples)
            ['sample_0' 'sample_1' 'sample_2' 'sample_3']

        """
        return self._metadata2_memmaps["samples"]

    @property
    def ids(self) -> List[str]:
        """
        The variant identifiers (a :class:`numpy.ndarray` of ``str``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.ids)
            ['SNP1' 'SNP2' 'SNP3' 'SNP4']

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     ids = bgen.ids.copy()
            >>> print(ids)
            ['SNP1' 'SNP2' 'SNP3' 'SNP4']


        """
        return self._metadata2_memmaps["ids"]

    @property
    def rsids(self) -> List[str]:
        """
        The variant RS numbers (a :class:`numpy.ndarray` of ``str``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.rsids)
            ['RS1' 'RS2' 'RS3' 'RS4']

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     rsids = bgen.rsids.copy()
            >>> print(rsids)
            ['RS1' 'RS2' 'RS3' 'RS4']

        """
        return self._metadata2_memmaps["rsids"]

    @property
    def _vaddr(self) -> List[int]:
        return self._metadata2_memmaps["vaddr"]

    @property
    def _ncombinations(self) -> List[int]:
        return self._metadata2_memmaps["ncombinations"]

    @property
    def _sample_range(self) -> List[int]:
        return self._metadata2_memmaps["sample_range"]

    @property
    def chromosomes(self) -> List[str]:
        """
        The chromosome of each variant (a :class:`numpy.ndarray` of ``str``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.chromosomes)
            ['1' '1' '1' '1']

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     chromosomes = bgen.chromosomes.copy()
            >>> print(chromosomes)
            ['1' '1' '1' '1']

        """
        return self._metadata2_memmaps["chromosomes"]

    @property
    def positions(self) -> List[int]:
        """
        The genetic position of each variant (a :class:`numpy.ndarray` of ``int``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.positions)
            [1 2 3 4]

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     positions = bgen.positions.copy()
            >>> print(positions)
            [1 2 3 4]

        """
        return self._metadata2_memmaps["positions"]

    @property
    def nalleles(self) -> List[int]:
        """
        The number of alleles for each variant (a :class:`numpy.ndarray` of ``int``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.nalleles)
            [2 2 2 2]

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     nalleles = bgen.nalleles.copy()
            >>> print(nalleles)
            [2 2 2 2]

        """
        return self._metadata2_memmaps["nalleles"]

    @property
    def allele_ids(self) -> List[str]:
        """
        The comma-delimited list of alleles for each variant (a :class:`numpy.ndarray` of ``str``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.allele_ids)
            ['A,G' 'A,G' 'A,G' 'A,G']

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     allele_ids = bgen.allele_ids.copy()
            >>> print(allele_ids)
            ['A,G' 'A,G' 'A,G' 'A,G']

        """
        return self._metadata2_memmaps["allele_ids"]

    @property
    def ncombinations(self) -> List[int]:
        """
        The number of values needed for each variant's probability distribution (a
        :class:`numpy.ndarray` of ``int``).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.ncombinations)
            [4 4 4 4]

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     ncombinations = bgen.ncombinations.copy()
            >>> print(ncombinations)
            [4 4 4 4]

        """
        return self._metadata2_memmaps["ncombinations"]

    @property
    def phased(self) -> List[bool]:
        """
        For each variant, ``True`` if and only the variant is phased (a :class:`numpy.ndarray` of
        bool).

        Note
        ----
        To access after file closes, make a copy.

        Example
        -------
        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> file = example_filepath("haplotypes.bgen")
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     print(bgen.phased)
            [ True  True  True  True]

            >>> # To access after file closes, make a copy
            >>> with open_bgen(file, verbose=False) as bgen:
            ...     phased = bgen.phased.copy()
            >>> print(phased)
            [ True  True  True  True]

        """
        return self._metadata2_memmaps["phased"]

    @staticmethod
    def _split_index(index):
        if not isinstance(index, tuple):
            index = (None, index)
        samples_index = open_bgen._fix_up_index(index[0])
        variants_index = open_bgen._fix_up_index(index[1])
        return samples_index, variants_index

    @staticmethod
    def _fix_up_index(index):
        if index is None:  # make a shortcut for None
            return slice(None)
        try:  # If index is an int, return it in an array
            index = index.__index__()  # (see
            # https://stackoverflow.com/questions/3501382/checking-whether-a-variable-is-an-integer-or-not)
            return [index]
        except Exception:
            pass
        return index

    def _extract_nalleles_ids_etc(self, mmm_wplus):
        with tmp_cwd():
            metafile_filepath = Path("bgen.metadata")
            self._cbgen.create_metafile(metafile_filepath, verbose=self._verbose)

            with _log_in_place("metadata", self._verbose) as updater:
                with bgen_metafile(metafile_filepath) as mf:
                    nparts = mf.npartitions

                    vaddr_memmap = mmm_wplus.append_empty(
                        "vaddr", (self.nvariants), "uint64"
                    )
                    positions_memmap = mmm_wplus.append_empty(
                        "positions", (self.nvariants), "uint32"
                    )
                    nalleles_memmap = mmm_wplus.append_empty(
                        "nalleles", (self.nvariants), "uint16"
                    )

                    vid_max_max = 0
                    rsid_max_max = 0
                    chrom_max_max = 0
                    allele_ids_max_max = 0

                    start = 0
                    for ipart2 in range(nparts):  # LATER multithread?
                        # LATER in notebook this message doesn't appear on one line
                        updater(
                            "'nallele': part {0:,} of {1:,}".format(ipart2 + 1, nparts)
                        )

                        partition = mf.read_partition(ipart2)
                        variants = partition.variants
                        nvariants = partition.variants.size

                        vid_max_max = max(vid_max_max, variants.id.dtype.itemsize)
                        rsid_max_max = max(rsid_max_max, variants.rsid.dtype.itemsize)
                        chrom_max_max = max(
                            chrom_max_max, variants.chromosome.dtype.itemsize
                        )
                        allele_ids_max_max = max(
                            allele_ids_max_max, variants.allele_ids.dtype.itemsize
                        )

                        end = start + nvariants
                        vaddr_memmap[start:end] = variants.offset
                        positions_memmap[start:end] = variants.position
                        nalleles_memmap[start:end] = variants.nalleles
                        start = end

                    ids_memmap = mmm_wplus.append_empty(
                        "ids", (self.nvariants), f"<U{vid_max_max}"
                    )
                    rsids_memmap = mmm_wplus.append_empty(
                        "rsids", (self.nvariants), f"<U{rsid_max_max}"
                    )
                    chrom_memmap = mmm_wplus.append_empty(
                        "chromosomes", (self.nvariants), f"<U{chrom_max_max}"
                    )
                    allele_ids_memmap = mmm_wplus.append_empty(
                        "allele_ids", (self.nvariants), f"<U{allele_ids_max_max}"
                    )

                    start = 0
                    for ipart2 in range(nparts):  # LATER multithread?
                        updater("'ids': part {0:,} of {1:,}".format(ipart2 + 1, nparts))

                        partition = mf.read_partition(ipart2)
                        variants = partition.variants
                        nvariants = partition.variants.size

                        end = start + nvariants
                        ids_memmap[start:end] = variants.id
                        rsids_memmap[start:end] = variants.rsid
                        chrom_memmap[start:end] = variants.chromosome
                        allele_ids_memmap[start:end] = variants.allele_ids
                        start = end

    def __str__(self):
        return "{0}('{1}')".format(self.__class__.__name__, self._filepath.name)

    def close(self):
        """
        Close a :class:`open_bgen` object that was opened for reading.

        Notes
        -----
        Better alternatives to :meth:`close` include the
        `with <https://docs.python.org/3/reference/compound_stmts.html#grammar-token-with-stmt>`__
        statement (closes the file automatically) and the `del
        <https://docs.python.org/3/reference/simple_stmts.html#grammar-token-del-stmt>`__
        statement (which closes the file and *deletes* the object).
        Doing nothing, while not better, is usually fine.

        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>> file = example_filepath("haplotypes.bgen")
            >>> bgen = open_bgen(file, verbose=False)
            >>> print(bgen.read(2))
            [[[1. 0. 0. 1.]]
            <BLANKLINE>
             [[0. 1. 0. 1.]]
            <BLANKLINE>
             [[1. 0. 1. 0.]]
            <BLANKLINE>
             [[0. 1. 1. 0.]]]
            >>> bgen.close()     #'del bgen' is better.

        """
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        if (
            hasattr(self, "_cbgen") and self._cbgen is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            self._cbgen.close()
            del self._cbgen
            # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect.
        if (
            hasattr(self, "_metadata2_memmaps") and self._metadata2_memmaps is not None
        ):  # we need to test this because Python doesn't guarantee that __init__ was
            # fully run
            self._metadata2_memmaps.__exit__(None, None, None)
            del (
                self._metadata2_memmaps
            )  # This allows __del__ and __exit__ to be called twice on the same object with
            # no bad effect.

    def allele_expectation(
        self,
        index: Optional[Any] = None,
        assume_constant_ploidy: bool = True,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Allele expectation.

        Parameters
        ----------
        index
            An expression specifying the samples and variants of interest. (See :ref:`read_examples` in :meth:`.read` for details.)
            Defaults to ``None``, meaning compute for all samples and variants.
        assume_constant_ploidy: bool
            When ploidy count can be assumed to be constant, calculations are much faster.
            Defaults to ``True``.


        Returns
        -------
        one or two :class:`numpy.ndarray`:
            always in this order

            * Samples-by-variants-by-alleles matrix of allele expectations,
            * Samples-by-variants-by-alleles matrix of frequencies, if ``return_frequencies`` is ``True``

        Note
        ----
        This method supports unphased genotypes only.


        .. _allele_expectation_examples:

        Examples
        --------
        .. doctest::

            >>> from bgen_reader import allele_expectation, example_filepath, read_bgen
            >>> from texttable import Texttable
            >>>
            >>> filepath = example_filepath("example.32bits.bgen")
            >>>
            >>> # Read the example.
            >>> bgen = open_bgen(filepath, verbose=False)
            >>> sample_index = bgen.samples=="sample_005" # will be only 1 sample
            >>> variant_index = bgen.rsids=="RSID_6"      # will be only 1 variant
            >>> p = bgen.read((sample_index,variant_index))
            >>> # Allele expectation makes sense for unphased genotypes only,
            >>> # which is the case here.
            >>> e = bgen.allele_expectation((sample_index,variant_index))
            >>> alleles_per_variant = [allele_ids.split(',') for allele_ids in bgen.allele_ids[variant_index]]
            >>>
            >>> # Print what we have got in a nice format.
            >>> table = Texttable()
            >>> table = table.add_rows(
            ...     [
            ...         ["", "AA", "AG", "GG", "E[.]"],
            ...         ["p"] + list(p[0,0,:]) + ["na"],
            ...         ["#" + alleles_per_variant[0][0], 2, 1, 0, e[0,0,0]],
            ...         ["#" + alleles_per_variant[0][1], 0, 1, 2, e[0,0,1]],
            ...     ]
            ... )
            >>> print(table.draw())
            +----+-------+-------+-------+-------+
            |    |  AA   |  AG   |  GG   | E[.]  |
            +====+=======+=======+=======+=======+
            | p  | 0.012 | 0.987 | 0.001 | na    |
            +----+-------+-------+-------+-------+
            | #A | 2     | 1     | 0     | 1.011 |
            +----+-------+-------+-------+-------+
            | #G | 0     | 1     | 2     | 0.989 |
            +----+-------+-------+-------+-------+


        If ``return_frequencies`` is true, this method will also return the allele frequency.

        .. doctest::

            >>> from bgen_reader import open_bgen, example_filepath
            >>>
            >>> filepath = example_filepath("example.32bits.bgen")
            >>> bgen = open_bgen(filepath, verbose=False)
            >>>
            >>> variant_index = (bgen.rsids=="RSID_6")      # will be only 1 variant
            >>> e = bgen.allele_expectation(variant_index)
            >>> f = bgen.allele_frequency(e)
            >>> alleles_per_variant = [allele_ids.split(',') for allele_ids in bgen.allele_ids[variant_index]]
            >>> print(alleles_per_variant[0][0] + ": {}".format(f[0,0]))
            A: 229.23103218810434
            >>> print(alleles_per_variant[0][1] + ": {}".format(f[0,1]))
            G: 270.7689678118956
            >>> print(bgen.ids[variant_index][0],bgen.rsids[variant_index][0])
            SNPID_6 RSID_6

        To find dosage, just select the column of interest from the expectation.

        .. doctest::

            >>> from bgen_reader import example_filepath, open_bgen
            >>>
            >>> filepath = example_filepath("example.32bits.bgen")
            >>>
            >>> # Read the example.
            >>> bgen = open_bgen(filepath, verbose=False)
            >>>
            >>> # Extract the allele expectations of the fourth variant.
            >>> variant_index = 3
            >>> e = bgen.allele_expectation(variant_index)
            >>>
            >>> # Compute the dosage when considering the allele
            >>> # in position 1 as the reference/alternative one.
            >>> alt_allele_index = 1
            >>> dosage = e[...,1]
            >>>
            >>> # Print the dosage for only the first five samples
            >>> # and the one (and only) variant
            >>> print(dosage[:5,0])
            [1.96185308 0.00982666 0.01745552 1.00347899 1.01153563]
            >>> del bgen
            >>>
            >>> import pandas as pd
            >>> from bgen_reader import open_bgen
            >>> filepath = example_filepath("example.32bits.bgen")
            >>> bgen = open_bgen(filepath, verbose=False)
            >>>
            >>> variant_index = [3]
            >>> # Print the metadata of the fourth variant.
            >>> print(bgen.ids[variant_index],bgen.rsids[variant_index])
            ['SNPID_5'] ['RSID_5']
            >>> probs, missing, ploidy = bgen.read(variant_index,return_missings=True,return_ploidies=True)
            >>> print(np.unique(missing),np.unique(ploidy))
            [False] [2]
            >>> df1 = pd.DataFrame({'sample':bgen.samples,'0':probs[:,0,0],'1':probs[:,0,1],'2':probs[:,0,2]})
            >>> print(df1) # doctest: +NORMALIZE_WHITESPACE
                        sample        0        1        2
            0    sample_001  0.00488  0.02838  0.96674
            1    sample_002  0.99045  0.00928  0.00027
            2    sample_003  0.98932  0.00391  0.00677
            3    sample_004  0.00662  0.98328  0.01010
            ..          ...      ...      ...      ...
            496  sample_497  0.00137  0.01312  0.98550
            497  sample_498  0.00552  0.99423  0.00024
            498  sample_499  0.01266  0.01154  0.97580
            499  sample_500  0.00021  0.98431  0.01547
            <BLANKLINE>
            [500 rows x 4 columns]
            >>> alleles_per_variant = [allele_ids.split(',') for allele_ids in bgen.allele_ids[variant_index]]
            >>> e = bgen.allele_expectation(variant_index)
            >>> f = bgen.allele_frequency(e)
            >>> df2 = pd.DataFrame({'sample':bgen.samples,alleles_per_variant[0][0]:e[:,0,0],alleles_per_variant[0][1]:e[:,0,1]})
            >>> print(df2)  # doctest: +NORMALIZE_WHITESPACE
                        sample        A        G
            0    sample_001  0.03815  1.96185
            1    sample_002  1.99017  0.00983
            2    sample_003  1.98254  0.01746
            3    sample_004  0.99652  1.00348
            ..          ...      ...      ...
            496  sample_497  0.01587  1.98413
            497  sample_498  1.00528  0.99472
            498  sample_499  0.03687  1.96313
            499  sample_500  0.98474  1.01526
            <BLANKLINE>
            [500 rows x 3 columns]
            >>> df3 = pd.DataFrame({'allele':alleles_per_variant[0],bgen.rsids[variant_index][0]:f[0,:]})
            >>> print(df3)
              allele    RSID_5
            0      A 305.97218
            1      G 194.02782
            >>> alt_index = f[0,:].argmin()
            >>> alt = alleles_per_variant[0][alt_index]
            >>> dosage = e[:,0,alt_index]
            >>> df4 = pd.DataFrame({'sample':bgen.samples,f"alt={alt}":dosage})
            >>> # Dosages when considering G as the alternative allele.
            >>> print(df4) # doctest: +NORMALIZE_WHITESPACE
                     sample    alt=G
            0    sample_001  1.96185
            1    sample_002  0.00983
            2    sample_003  0.01746
            3    sample_004  1.00348
            ..          ...      ...
            496  sample_497  1.98413
            497  sample_498  0.99472
            498  sample_499  1.96313
            499  sample_500  1.01526
            <BLANKLINE>
            [500 rows x 2 columns]

        """
        _, variants_index = self._split_index(index)
        phased_list = self.phased[variants_index]
        nalleles = self.nalleles[variants_index]
        if any(phased_list):
            raise ValueError(
                "Allele expectation is define for unphased genotypes only."
            )

        if len(np.unique(nalleles)) > 1:
            raise ValueError(
                "Current code requires that all selected variants have the same number of alleles"
            )
        if assume_constant_ploidy:
            ploidy0 = self.read(return_probabilities=False, return_ploidies=True)[
                [0], 0
            ]
            genotype = get_genotypes(ploidy0, self.nalleles[0])[0]
            count = asarray(genotypes_to_allele_counts(genotype), float)

            probs = self.read(index)
            if (
                np.product(probs.shape[:2]) == 0
            ):  # handle the case where user asks for no samples or no variants
                expecx = np.zeros((probs.shape[0], probs.shape[1], count.shape[-1]))
            else:
                if probs.shape[-1] != count.shape[0]:
                    raise ValueError("Try 'assume_constant_ploidy=False'")
                expecx = probs.dot(count)
        else:
            probs, ploidy = self.read(index, return_ploidies=True)
            if (
                np.product(probs.shape[:2]) == 0
            ):  # handle the case where user asks for no samples or no variants
                expecx = np.zeros((probs.shape[0], probs.shape[1], self.nalleles[0]))
            else:
                outer_expec = []
                for vi in range(probs.shape[1]):  # for each variant ...
                    genotypes = get_genotypes(ploidy[:, vi], nalleles[vi])
                    probsvi = probs[:, vi, :]
                    expec = []
                    for i, genotype in enumerate(genotypes):
                        count = asarray(genotypes_to_allele_counts(genotype), float)
                        n = count.shape[0]
                        expec.append((count.T * probsvi[i, :n]).sum(1))
                    outer_expec.append(stack(expec, axis=0))
                expecx = stack(outer_expec, axis=1)

        return expecx

    @staticmethod
    def allele_frequency(allele_expectation: np.ndarray) -> np.ndarray:
        """
        Allele expectation frequency.

        You have to provide the allele expectations, :meth:`.allele_expectation`.
        """
        nallele0 = allele_expectation.shape[-1]
        return allele_expectation.sum(0) / nallele0

    def __del__(self):
        self.__exit__()
