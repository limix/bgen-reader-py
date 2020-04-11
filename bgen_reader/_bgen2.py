import os
import numpy as np
from contextlib import contextmanager
from tempfile import mkdtemp
import shutil
from pathlib import Path
from typing import Optional, Union
from numpy import empty, uint16, uint32, uint64, zeros

#!!!cmk figure out a good name of this file
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

#!!!cmk add type info
#!!!cmk write doc
#!!!cmk test doc
#!!!cmk1 create test cases that generate big datasets (if qctool is available)
#!!!cmk1 create test cases with full coverage
#!!!cmk ok to have metadata2.npz location be fixed for now?
class open_bgen(object):
    def __init__(self, filepath: Union[str, Path],
                 samples_filepath: Optional[Union[str, Path]] = None,
                 verbose : bool = False):
        filepath = Path(filepath)
        assert_file_exist(filepath)
        assert_file_readable(filepath)


        self._verbose = verbose
        self._filepath = filepath

        self._bgen_context_manager = bgen_file(filepath)
        self._bgen = self._bgen_context_manager.__enter__()

        self._samples = np.array(self._get_samples(samples_filepath),dtype='str')
        self._sample_range = np.arange(len(self._samples),dtype=np.int)

        metadata2 = filepath.with_suffix('.metadata2.npz')
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

    def _get_samples(self,sample_file): #!!!cmk similar code in _reader.py
        if sample_file is None:
            if self._bgen.contain_samples:
                return self._bgen.read_samples()
            else:
                return generate_samples(self._bgen.nsamples)
        else:
            samples_filepath = Path(sample_file)
            assert_file_exist(sample_file)
            assert_file_readable(sample_file)
            return read_samples_file(sample_file, self._verbose)

    @staticmethod
    def _fix_up_index(index):
        if type(index) is np.int: #!!!make sure this works with all int types
            return [index]
        return index


    #!!!cmk test each option
    def read(self, index=None, dtype=np.float, order='F', max_combinations=None, return_probabilities=True, return_missings=False, return_ploidies=False):
        '''
        !!!cmkwrite doc
        !!!cmk tell probs will be 3D array
        !!!cmk if both missing and ploidy are returned, missing will be first. and both will be 2-D arrays
        '''
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
        #!!!cmk if verbose is true, give some status
        for out_index,vaddr0 in enumerate(vaddr):
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
        with bgen_metafile(Path(metafile_filepath)) as mf:
            nparts = mf.npartitions
            id_list, rsid_list,chrom_list,position_list,vaddr_list,nalleles_list,allele_ids_list,ncombinations_list,phased_list = [],[],[],[],[],[],[],[],[]

            #!!!If verbose, should tell how it is going
            for ipart in range(nparts): #LATER multithread?

                #!!!cmk this code is very similar to other code
                partition = lib.bgen_metafile_read_partition(mf._bgen_metafile, ipart)
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
            if self._verbose and len(id_list)%1000==0:
                print('{0}'.format(len(id_list)))#!!!cmk put in standard-style status message. Elsewhere, too?
            genotype = lib.bgen_file_open_genotype(self._bgen._bgen_file, vaddr0)
            ncombinations_list.append(lib.bgen_genotype_ncombs(genotype))
            phased_list.append(lib.bgen_genotype_phased(genotype))
            lib.bgen_genotype_close(genotype)

        self._ncombinations = np.array(ncombinations_list,dtype='int')
        self._phased = np.array(phased_list,dtype='bool')

    @property
    def nvariants(self) -> int:
        '''!!!cmk doc this all all others
        '''
        return len(self.id)

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

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self._filepath)

    #!!!cmk here and elsewhere tell users to use 'del bgen12' not 'bgen12' so that object will be really gone and not just a zombie
    
    def __enter__(self):
        return self

    def __exit__(self, *_):
        if hasattr(self,'_bgen_context_manager') and self._bgen_context_manager is not None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self._bgen_context_manager.__exit__(None,None,None)
            del self._bgen_context_manager # This allows __del__ and __exit__ to be called twice on the same object with no bad effect.

    def __del__(self):
        self.__exit__()

#!!!cnj check how this works (if at all) with the other parts of the API (Dosage, Expectation, ...)   

if __name__ == "__main__":
    if True:
        #filepath = r'm:\deldir\1000x500000.bgen'
        #filepath = r'D:\OneDrive\Shares\bgenraaderpy\1x1000000.bgen'
        filepath = r'M:\del35\temp1024x16384-8.bgen'
        with open_bgen(filepath) as bgen2:
            print(bgen2.ids[:5]) #other property arrays include risd,chromosomes,positions,nallels, and allele_ids
            geno = bgen2.read(-1) # read the 200,000th variate's data
            #geno = bgen2.read() # read all, uses the ncombinations from the first variant
            geno = bgen2.read(slice(5)) # read first 5, uses the ncombinations from the first variant
            geno = bgen2.read(bgen2.chromosomes=='5',max_combinations=4) # read chrom1, set max_combs explicitly
            geno, missing = bgen2.read(0,return_missings=True)
        bgen2 = open_bgen(filepath)
        geno, missing = bgen2.read(0,return_missings=True)
    if True:
        filepath = r'm:\deldir\2500x500000.bgen'
        with open_bgen(filepath,verbose=True) as bgen2:
            print(bgen2.read(0)[0,0,:])
            print(bgen2.read(-1)[0,0,:])
    print('!!!done')
