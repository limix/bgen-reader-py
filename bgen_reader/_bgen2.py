import os
import numpy as np
from contextlib import contextmanager
from tempfile import mkdtemp
import shutil
from pathlib import Path
from typing import Optional, Union
from numpy import empty, uint16, uint32, uint64, zeros
import time
import datetime
import sys

#!!!cmk change 'from bgen_reader.' to 'from .'
from bgen_reader._bgen_file import bgen_file
from bgen_reader._samples import generate_samples, read_samples_file
from bgen_reader._file import (
    assert_file_exist,
    assert_file_readable,
    is_file_writable,
    path_to_filename,
)
from bgen_reader._bgen_metafile import bgen_metafile
from bgen_reader._ffi import ffi, lib
from bgen_reader._helper import _log_in_place
from bgen_reader.test.write_random import _write_random

#!!!cmk add type info
#!!!cmk write doc
#!!!cmk test doc
#!!!cmk ok to have metadata2.npz location be fixed for now?
class open_bgen(object):
    """
    Open a BGEN file for reading.

    Parameters
    ----------
    filepath
        Bgen file path.
    samples_filepath
        Path to a `sample format`_ file or ``None`` to read samples from the bgen file itself.
        Defaults to ``None``.
    verbose
        ``True`` to show progress; ``False`` otherwise. Defaults to ``True``.

    Returns
    -------
    an :class:`open_bgen` object : :class:`open_bgen`

    Examples
    --------
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
        >>>
        >>> bgen = open_bgen(file, verbose=False)
        >>> print(bgen.read((0,1)))  # read 1st sample and 2nd variant
        [[[0. 1. 1. 0.]]]
        >>> del bgen                 # close and delete object

    .. _sample format: https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html
    """
    def __init__(self, filepath: Union[str, Path],
                 samples_filepath: Optional[Union[str, Path]] = None,
                 verbose : bool = True):
        filepath = Path(filepath)
        assert_file_exist(filepath)
        assert_file_readable(filepath)

        self._verbose = verbose
        self._filepath = filepath


        self._bgen_context_manager = bgen_file(filepath)
        self._bgen = self._bgen_context_manager.__enter__()

        self._samples = np.array(self._get_samples(samples_filepath),dtype='str')
        self._sample_range = np.arange(len(self._samples),dtype=np.int)

        metadata2 = self._metadatapath_from_filename(filepath)
        if metadata2.exists() and os.path.getmtime(metadata2) < os.path.getmtime(filepath):
            metadata2.unlink()
        if metadata2.exists():
            d = np.load(str(metadata2))
            self._ids = d['ids']
            self._rsids = d['rsids']
            self._vaddr = d['vaddr']
            self._chromosomes = d['chromosomes']
            self._positions = d['positions']
            self._nalleles = d['nalleles']
            self._allele_ids = d['allele_ids']
            self._ncombinations = d['ncombinations']
            self._phased = d['phased']
        else:
            tempdir = None
            try:
                tempdir = mkdtemp(prefix='pysnptools')
                metafile_filepath = tempdir+'/bgen.metadata'
                self._bgen.create_metafile(Path(metafile_filepath),verbose=self._verbose)
                self._map_metadata(metafile_filepath)
                np.savez(metadata2,ids=self._ids,rsids=self._rsids,vaddr=self._vaddr,chromosomes=self._chromosomes,positions=self._positions,
                         nalleles=self._nalleles,allele_ids=self._allele_ids,ncombinations=self._ncombinations,phased=self._phased)
            finally:
                if tempdir is not None:
                    shutil.rmtree(tempdir)

        self.max_combinations = max(self._ncombinations)

    #This is static so that test code can use it easily.        
    @staticmethod
    def _metadatapath_from_filename(filename):
        return Path(filename).with_suffix('.metadata2.npz')

    @property
    def samples(self):
        return self._samples

    @property
    def ids(self):
        return self._ids

    @property
    def rsids(self):
        return self._rsids

    @property
    def chromosomes(self):
        return self._chromosomes

    @property
    def positions(self):
        return self._positions

    @property
    def nalleles(self):
        return self._nalleles

    @property
    def allele_ids(self):
        return self._allele_ids

    #!!!cmk these need types and docstrings
    @property
    def ncombinations(self):
        return self._ncombinations

    @property
    def phased(self):
        return self._phased

    @property
    def nvariants(self) -> int:
        '''!!!cmk doc this all all others
        '''
        return len(self.ids)

    #!!!cmk add a close method and test that any readers afterword raise a ValueError

    @property
    def nsamples(self) -> int:
        '''!!!cmk doc this all all others
        '''
        return len(self._samples)

    @property
    def shape(self) -> (int,int,int):
        '''!!!cmk doc this all all others
        '''
        return (self.nsamples,self.nvariants,self.max_combinations)

    def _get_samples(self,sample_file): #!!!cmk similar code in _reader.py
        if sample_file is None:
            if self._bgen.contain_samples:
                return self._bgen.read_samples()
            else:
                return generate_samples(self._bgen.nsamples)
        else:
            samples_filepath = Path(sample_file)
            assert_file_exist(samples_filepath)
            assert_file_readable(samples_filepath)
            return read_samples_file(samples_filepath, self._verbose)

    @staticmethod
    def _fix_up_index(index):
        if index is None: #make a shortcut for None
            return index
        try: #If index is an int, return it in an array
            index = index.__index__() #(see https://stackoverflow.com/questions/3501382/checking-whether-a-variable-is-an-integer-or-not)
            return [index] #!!!cmk should we return one less dimension?
        except:
            pass
        return index

    def read(self, index=None,
                   dtype = np.float64, #!!!cmk : Union(type,str) can't get to work
                   order: str ='F',
                   max_combinations: int = None, #cmk Optional(int)
                   return_probabilities: bool = True,
                   return_missings: bool = False,
                   return_ploidies: bool = False):
        #!!!cmkwrite doc
        #!!!cmk tell probs will be 3D array
        #!!!cmk if both missing and ploidy are returned, missing will be first. and both will be 2-D arrays
        """
        Read from an open_bgen object.

        Parameters
        ----------
        index
            An expression specifying the samples and variants to read (see below). 
            Defaults to ``None``, meaning read all.
        dtype : data-type
            The desired data-type for the returned probability array. #cmk somewhere give avoid about bits and this
            Defaults to ``np.float64``.
        order : {'F','C'}
            The desired memory layout for the returned probability array.
            Defaults to ``F`` (Fortran order, variant-major)
        max_combinations : int
            The number of values to allocate for each probability distribution. For example, 3 for unphased diploid data and 4 for phased diploid data.
            Any overallocated space is filled with ``np.nan``.
            Defaults to the number just large enough for any data in the file.
        return_probabilities: bool
            Read and return the probabilities for samples and variants specified.
            Defaults to ``True``
        return_missings: bool
            Return a boolean array telling which probabilities are missing.
            Defaults to ``False``
        return_ploidies: bool
            Read and return the ploidy for the samples and variants specified.
            Defaults to ``False``

        Returns
        -------
        zero to three :class:`numpy.ndarray`
            if return_probabilities is ``True`` (the default), the first return value will be 
                a ``dtype`` array of size (nsamples_out,nvariants_out,max_combinations).
            if return_missings is ``True``, the next return value will be a bool array of size (nsamples_out,nvariants_out).
            if return_ploidies is ``True``, the next return value will be an int array of size (nsamples_out,nvariants_out).

        Index
        -----

        Read all values:

        ``None``: Read all samples and variants

        Read selected variants:
        *variant_index*          : Real all samples from the variant specified where *variant_index* is of the form:

        *ivariant* (an ``int``) : Read all samples from the *ivariant* variant.
        [*ivariant0*, ...,*ivarianty*] (an ``list of ints``) : Read all samples from the *ivariant*s given
        slice(*variant_start*,*variant_stop*,*variant_step*) : Read all samples from the slice of variants given.
        [*bool0*, ...,*booln*] (an ``list of bools``) : Read all samples from the variants where the ``bool`` is ``True``.

        Read selected samples:

        (*sample_index*,None) : Read all the variants for the sample's specified, where *sample_index* follows the form of *variant_index*.

        Read selected samples and variants:

        (*sample_index*,*variant_index*)


        Examples
        --------
        .. doctest::

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
            >>>
            >>> bgen_e = open_bgen(example_filepath("example.bgen"), verbose=False)
            >>> print(bgen_e.read(5)) #read all samples for variant at position 5
            [[[           nan            nan            nan]]
            <BLANKLINE>
                [[3.63170030e-03 9.94598226e-01 1.77007378e-03]]
            <BLANKLINE>
                [[4.54711821e-03 9.90600594e-01 4.85228794e-03]]
            <BLANKLINE>
                ...
            <BLANKLINE>
                [[1.37329078e-03 2.03552104e-02 9.78271499e-01]]
            <BLANKLINE>
                [[2.25524965e-02 9.77172846e-01 2.74657970e-04]]
            <BLANKLINE>
                [[6.68334728e-03 1.22069847e-04 9.93194583e-01]]]
            >>> print(bgen_e.read(-1)) #read all samples for the last variant
            [[[1.33971970e-02 9.81353784e-01 5.24901878e-03]]
            <BLANKLINE>
             [[2.53295994e-02 7.32422108e-04 9.73937978e-01]]
            <BLANKLINE>
             [[9.97924526e-03 9.82696538e-01 7.32421666e-03]]
            <BLANKLINE>
             ...
            <BLANKLINE>
             [[9.79919470e-01 1.10168052e-02 9.06372443e-03]]
            <BLANKLINE>
             [[3.32641951e-03 1.50756983e-02 9.81597882e-01]]
            <BLANKLINE>
             [[5.88379060e-02 9.30725093e-01 1.04370010e-02]]]
            >>> print(bgen_e.read([4,5,2])) #read all samples for variants at position 4,5 and 2.
            [[[2.92977947e-03 9.96612442e-01 4.57778107e-04]
              [           nan            nan            nan]
              [9.92095944e-01 3.05176014e-04 7.59887952e-03]]
            <BLANKLINE>
             [[9.93530268e-01 2.38037063e-03 4.08936106e-03]
              [3.63170030e-03 9.94598226e-01 1.77007378e-03]
              [1.22986045e-02 9.81567363e-01 6.13403227e-03]]
            <BLANKLINE>
             [[9.84344482e-01 7.14111352e-03 8.51440407e-03]
              [4.54711821e-03 9.90600594e-01 4.85228794e-03]
              [4.79125930e-03 1.10778981e-02 9.84130843e-01]]
            <BLANKLINE>
             ...
            <BLANKLINE>
             [[1.31226042e-02 1.22070080e-04 9.86755326e-01]
              [1.37329078e-03 2.03552104e-02 9.78271499e-01]
              [6.40869071e-03 3.66210938e-04 9.93225098e-01]]
            <BLANKLINE>
             [[2.07519927e-03 3.05175781e-04 9.97619625e-01]
              [2.25524965e-02 9.77172846e-01 2.74657970e-04]
              [3.90624884e-03 1.77001953e-03 9.94323732e-01]]
            <BLANKLINE>
             [[5.79833985e-04 3.54004046e-03 9.95880126e-01]
              [6.68334728e-03 1.22069847e-04 9.93194583e-01]
              [1.22070080e-04 3.72313918e-03 9.96154791e-01]]]
            >>>
            >>> print(bgen_e.read(slice(5)) #read all samples for first 5 variants
            [[[0,3,3]]]
            >>> print(bgen_e.read(slice(2,5)) #read all samples for variants from position 2 (inclusive) to 5 (exclusive)
            [[[0,3,3]]]
            >>> print(bgen_e.read(slice(2,None)) #read all samples for variants starting at position 2
            [[[0,3,3]]]
            >>> print(bgen_e.read(slice(None,None,10)) #read all samples for every 10th variant
            [[[0,3,3]]]
            >>> print(bgen_e.read(bgen_e.chromosomes=='5') #read all samples for all variants in chromosome 5
            [[[0,3,3]]]
            >>> print(bgen_e.read((0,None)) #read the first sample across all variants
            [[[0,3,3]]]
            >>> print(bgen_e.read((slice(10,20),slice(10)) #read for samples 10 (inclusive) to 20 (exclusive), read the first 10 variants.
            [[[0,3,3]]]

            >>> probs,missing,ploidy = bgen_e.read(return_missings=True,return_ploidies=True) #read probabilities, missingness, and ploidy
            >>> print(polidy)
        """
        #!!!cmk allow single ints, lists of ints, lists of bools, None, and slices
        #!!!cmk (DECIDE LATER) could allow strings (variant names) and lists of strings

        max_combinations = max_combinations or self.max_combinations #!!!cmk test user setting max_combinations to 0 and 1

        if not isinstance(index,tuple): #!!!test with np.s_[,]
            index = (None,index)
        #!!!cmk raise error if not 
        samples_index = self._fix_up_index(index[0])
        variants_index = self._fix_up_index(index[1])

        if samples_index is None:
            samples_index = self._sample_range
        else:
            samples_index = self._sample_range[samples_index]

        if variants_index is None:
            vaddr = self._vaddr
            ncombinations = self._ncombinations
        else:
            vaddr = self._vaddr[variants_index]
            ncombinations = self._ncombinations[variants_index]

        #allocating prob_buffer only when its size changes makes reading 10x5M data 30% faster
        if return_probabilities:
            val = np.full((len(samples_index), len(vaddr), max_combinations), np.nan, dtype=dtype, order=order) #!!!cmk test on selecting zero variants
            prob_buffer = None
        if return_missings:
            missing_val = np.full((len(samples_index), len(vaddr)), False, dtype='bool', order=order)
        if return_ploidies:
            ploidy_val = np.full((len(samples_index), len(vaddr)), 0, dtype='int', order=order)


        #LATER multithread?
        with _log_in_place('reading', self._verbose) as updater:
            for out_index,vaddr0 in enumerate(vaddr):
                if out_index%100==0: #!!!cmk improve the freq
                    updater('part {0:,} of {1:,}'.format(out_index,len(vaddr))) #!!!cmk make this nice

                genotype = lib.bgen_file_open_genotype(self._bgen._bgen_file, vaddr0)

                if return_probabilities:
                    if prob_buffer is None or ncombinations[out_index] != prob_buffer.shape[-1]:
                        prob_buffer = np.full((len(self._samples), ncombinations[out_index]), np.nan, order='C', dtype='float64')
                    lib.bgen_genotype_read(genotype, ffi.cast("double *", prob_buffer.ctypes.data))
                    val[:,out_index,:ncombinations[out_index]] = prob_buffer if (samples_index is self._sample_range) else prob_buffer[samples_index,:]

                if return_missings:
                    missing_val[:,out_index] = [lib.bgen_genotype_missing(genotype, i) for i in samples_index]

                if return_ploidies:
                    ploidy_val[:,out_index] = [lib.bgen_genotype_ploidy(genotype, i) for i in samples_index]

                lib.bgen_genotype_close(genotype)

        result_array = ([val] if return_probabilities else []) + ([missing_val] if return_missings else []) + ([ploidy_val] if return_ploidies else [])
        if len(result_array)==1:
            return result_array[0]
        else:
            return tuple(result_array)

    def _map_metadata(self,metafile_filepath): 
        with _log_in_place('metadata', self._verbose) as updater:
            with bgen_metafile(Path(metafile_filepath)) as mf:
                nparts = mf.npartitions
                id_list, rsid_list,chrom_list,position_list,vaddr_list,nalleles_list,allele_ids_list,ncombinations_list,phased_list = [],[],[],[],[],[],[],[],[]

                #!!!If verbose, should tell how it is going
                for ipart2 in range(nparts): #LATER multithread?
                    updater('step 2: part {0:,} of {1:,}'.format(ipart2,nparts)) #!!!cmk make this nice

                    #!!!cmk this code is very similar to other code
                    partition = lib.bgen_metafile_read_partition(mf._bgen_metafile, ipart2)
                    if partition == ffi.NULL:
                        raise RuntimeError(f"Could not read partition {partition}.")

                    #cmk similar code in _bgen_metafile.py
                    nvariants = lib.bgen_partition_nvariants(partition)

                    position = empty(nvariants, dtype=uint32)
                    nalleles = empty(nvariants, dtype=uint16)
                    offset = empty(nvariants, dtype=uint64)
                    vid_max_len = ffi.new("uint32_t[]", 1)
                    rsid_max_len = ffi.new("uint32_t[]", 1)
                    chrom_max_len = ffi.new("uint32_t[]", 1)
                    allele_ids_max_len = ffi.new("uint32_t[]", 1)
                    position_ptr = ffi.cast("uint32_t *", ffi.from_buffer(position))
                    nalleles_ptr = ffi.cast("uint16_t *", ffi.from_buffer(nalleles))
                    offset_ptr = ffi.cast("uint64_t *", ffi.from_buffer(offset))
                    lib.read_partition_part1(
                        partition,
                        position_ptr,
                        nalleles_ptr,
                        offset_ptr,
                        vid_max_len,
                        rsid_max_len,
                        chrom_max_len,
                        allele_ids_max_len,
                    )
                    vid = zeros(nvariants, dtype=f"S{vid_max_len[0]}")
                    rsid = zeros(nvariants, dtype=f"S{rsid_max_len[0]}")
                    chrom = zeros(nvariants, dtype=f"S{chrom_max_len[0]}")
                    allele_ids = zeros(nvariants, dtype=f"S{allele_ids_max_len[0]}")
                    lib.read_partition_part2(
                        partition,
                        ffi.from_buffer("char[]", vid),
                        vid_max_len[0],
                        ffi.from_buffer("char[]", rsid),
                        rsid_max_len[0],
                        ffi.from_buffer("char[]", chrom),
                        chrom_max_len[0],
                        ffi.from_buffer("char[]", allele_ids),
                        allele_ids_max_len[0],
                     )

                    id_list.append(vid)
                    rsid_list.append(rsid)
                    chrom_list.append(chrom)
                    position_list.append(position)
                    nalleles_list.append(nalleles)
                    allele_ids_list.append(allele_ids)
                    vaddr_list.append(offset)

            ##!!cmk do these need to pre-declared.                    
            #!!!cmk use concatenate(...out=) instead
            self._ids = np.array(np.concatenate(id_list),dtype='str') # dtype needed to make unicode
            self._rsids = np.array(np.concatenate(rsid_list),dtype='str')
            self._vaddr = np.concatenate(vaddr_list)
            self._chromosomes = np.array(np.concatenate(chrom_list),dtype='str') 
            self._positions = np.concatenate(position_list)
            self._nalleles = np.concatenate(nalleles_list)
            self._allele_ids = np.array(np.concatenate(allele_ids_list),dtype='str') #cmk check that main api doesn't return bytes

            for i,vaddr0 in enumerate(self._vaddr):
                if i%1000==0: #!!!cmk improve the freq
                    updater('step 3: part {0:,} of {1:,}'.format(i,self.nvariants)) #!!!cmk make this nice
                genotype = lib.bgen_file_open_genotype(self._bgen._bgen_file, vaddr0)
                ncombinations_list.append(lib.bgen_genotype_ncombs(genotype))
                phased_list.append(lib.bgen_genotype_phased(genotype))
                lib.bgen_genotype_close(genotype)

            self._ncombinations = np.array(ncombinations_list,dtype='int')
            self._phased = np.array(phased_list,dtype='bool')

    def __str__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self._filepath.name)

    #!!!cmk here and elsewhere tell users to use 'del bgen12' not 'bgen12' so that object will be really gone and not just a zombie
    
    def __enter__(self):
        return self

    def __exit__(self, *_):
        if hasattr(self,'_bgen_context_manager') and self._bgen_context_manager is not None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self._bgen_context_manager.__exit__(None,None,None)
            del self._bgen_context_manager # This allows __del__ and __exit__ to be called twice on the same object with no bad effect.

    def __del__(self):
        self.__exit__()

#!!!cmk check how this works (if at all) with the other parts of the API (Dosage, Expectation, ...)   



if __name__ == "__main__":
    if False:
        from bgen_reader.test.test_bgen2 import test_bigfile
        test_bigfile(verbose=True)
    if False:
        filepath = r'm:\deldir\400x500x3.bgen'
        _write_random(filepath, 400, 500, chrom_count=5, cleanup_temp_files=False)
    if False:
        #filepath = r'm:\deldir\1000x500000.bgen'
        #filepath = r'D:\OneDrive\Shares\bgenraaderpy\1x1000000.bgen'
        filepath = r'M:\del35\temp1024x16384-8.bgen'
        with open_bgen(filepath,verbose=True) as bgen2:
            print(bgen2.ids[:5]) #other property arrays include risd,chromosomes,positions,nallels, and allele_ids
            geno = bgen2.read(-1) # read the 200,000th variate's data
            #geno = bgen2.read() # read all, uses the ncombinations from the first variant
            geno = bgen2.read(slice(5)) # read first 5, uses the ncombinations from the first variant
            geno = bgen2.read(bgen2.chromosomes=='5',max_combinations=4) # read chrom1, set max_combs explicitly
            geno, missing = bgen2.read(0,return_missings=True)
        bgen2 = open_bgen(filepath)
        geno, missing = bgen2.read(0,return_missings=True)
    if False:
        filepath = r'm:\deldir\2500x500000.bgen'
        with open_bgen(filepath,verbose=True) as bgen2:
            print(bgen2.read(0)[0,0,:])
            print(bgen2.read(-1)[0,0,:])

    if True:
        import pytest

        pytest.main([__file__])
    print('!!!cmk')
