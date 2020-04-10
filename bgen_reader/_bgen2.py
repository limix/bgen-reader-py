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

#!!!cmk1 can we have it both a contextmanger and not?
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
        self.filepath = filepath

        self._bgen_context_manager = bgen_file(filepath)
        self._bgen = self._bgen_context_manager.__enter__()

        self.samples = np.array(self._get_samples(samples_filepath),dtype='str')

        metadata2 = filepath.with_suffix('.metadata2.npz')
        if metadata2.exists():
            d = np.load(str(metadata2))
            #!!!cmk why are some properties plural and others aren't?
            self.id = d['id']
            self.rsid = d['rsid']
            self.vaddr = d['vaddr']
            self.chrom = d['chrom']
            self.position = d['position']
            self.nalleles = d['nalleles']
            self.allele_ids = d['allele_ids']
            self.ncombs = d['ncombs']
            self.phased = d['phased']
        else:
            tempdir = None
            try:
                tempdir = mkdtemp(prefix='pysnptools')
                metafile_filepath = tempdir+'/bgen.metadata'
                self._bgen.create_metafile(Path(metafile_filepath),verbose=self._verbose)
                self._map_metadata(metafile_filepath)
                np.savez(metadata2,id=self.id,rsid=self.rsid,vaddr=self.vaddr,chrom=self.chrom,position=self.position,
                         nalleles=self.nalleles,allele_ids=self.allele_ids,ncombs=self.ncombs,phased=self.phased)
            finally:
                if tempdir is not None:
                    shutil.rmtree(tempdir)

        self.max_ncombs = max(self.ncombs)

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

    #!!!cmk test each option
    def read(self, variants=None, max_ncombs=None, dtype=np.float, order='F', return_missing=False, return_ploidy=False): #!!!cmk also allow samples to be selected?
        '''
        !!!cmkwrite doc
        !!!cmk tell probs will be 3D array
        !!!cmk if both missing and ploidy are returned, missing will be first. and both will be 2-D arrays
        '''
        #!!!cmk allow single ints, lists of ints, lists of bools, None, and slices
        #!!!cmk (DECIDE LATER) could allow strings (variant names) and lists of strings

        max_ncombs = max_ncombs or self.max_ncombs #!!!cmk test user setting max_ncombs to 0 and 1

        if type(variants) is np.int: #!!!make sure this works with all int types
            variants = [variants]
        if variants is None:
            vaddr = self.vaddr
            ncombs = self.ncombs
        else:
            vaddr = self.vaddr[variants]
            ncombs = self.ncombs[variants]

        #allocating p only once make reading 10x5M data 30% faster
        val = np.full((len(self.samples), len(vaddr), max_ncombs), np.nan, dtype=dtype, order=order) #!!!cmk test on selecting zero variants
        if return_missing:
            missing_val = np.full((len(self.samples), len(vaddr)), False, dtype='bool', order=order)
        if return_ploidy:
            ploidy_val = np.full((len(self.samples), len(vaddr)), 0, dtype='int', order=order)

        p = None #!!!cmk 'p' seems hard to read about 'prop_buffer'?

        #LATER multithread?
        #!!!cmk if verbose is true, give some status
        for out_index,vaddr0 in enumerate(vaddr):
            if p is None or ncombs[out_index] != p.shape[-1]:
                p = np.full((len(self.samples), ncombs[out_index]), np.nan, order='C', dtype='float64')
            genotype = lib.bgen_file_open_genotype(self._bgen._bgen_file, vaddr0)
            lib.bgen_genotype_read(genotype, ffi.cast("double *", p.ctypes.data))
            val[:,out_index,:ncombs[out_index]] = p

            if return_missing:
                missing_val[:,out_index] = [lib.bgen_genotype_missing(genotype, i) for i in range(self.nsamples)]
            if return_ploidy:
                ploidy_val[:,out_index] = [lib.bgen_genotype_ploidy(genotype, i) for i in range(self.nsamples)]

            lib.bgen_genotype_close(genotype)

        if not return_missing and not return_ploidy:
            return val
        else:
            val_array = [val] + ([missing_val] if return_missing else []) + ([ploidy_val] if return_ploidy else [])
            return tuple(val_array)

    def _map_metadata(self,metafile_filepath): 
        with bgen_metafile(Path(metafile_filepath)) as mf:
            nparts = mf.npartitions
            id_list, rsid_list,chrom_list,position_list,vaddr_list,nalleles_list,allele_ids_list,ncombs_list,phased_list = [],[],[],[],[],[],[],[],[]

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

        #!!!cmk do these need to be @properties?
        ##!!cmk do these need to pre-declared.                    
        #!!!cmk use concatenate(...out=) instead
        self.id = np.array(np.concatenate(id_list),dtype='str')#dtype needed to make unicode
        self.rsid = np.array(np.concatenate(rsid_list),dtype='str')
        self.vaddr = np.concatenate(vaddr_list)
        self.chrom = np.array(np.concatenate(chrom_list),dtype='str') 
        self.position = np.concatenate(position_list)
        self.nalleles = np.concatenate(nalleles_list)
        self.allele_ids = np.array(np.concatenate(nalleles_list),dtype='str') #cmk check that main api doesn't return bytes

        for i,vaddr0 in enumerate(self.vaddr):
            if self._verbose and len(id_list)%1000==0:
                print('{0}'.format(len(id_list)))#!!!cmk put in standard-style status message. Elsewhere, too?
            genotype = lib.bgen_file_open_genotype(self._bgen._bgen_file, vaddr0)
            ncombs_list.append(lib.bgen_genotype_ncombs(genotype))
            phased_list.append(lib.bgen_genotype_phased(genotype))
            lib.bgen_genotype_close(genotype)

        self.ncombs = np.array(ncombs_list,dtype='int')
        self.phased = np.array(phased_list,dtype='bool')

    @property
    def nvariants(self) -> int:
        '''!!!cmk doc this all all others
        '''
        return len(self.id)

    @property
    def nsamples(self) -> int:
        '''!!!cmk doc this all all others
        '''
        return len(self.samples)

    @property
    def shape(self) -> (int,int,int):
        '''!!!cmk doc this all all others
        '''
        return (self.nsamples,self.nvariants,self.max_ncombs)

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self.filepath)

    #!!!cmk here and elsewhere tell users to use 'del bgen12' not 'bgen12' so that object will be really gone and not just a zombie
    
    def __enter__(self):
        return self

    def __exit__(self, *_):
        if hasattr(self,'_bgen_context_manager') and self._bgen_context_manager is not None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self._bgen_context_manager.__exit__(None,None,None)
            del self._bgen_context_manager # This allows __del__ and __exit__ to be called twice on the same object with no bad effect.

    def __del__(self):
        self.__exit__()

if __name__ == "__main__":
    if True:
        #filepath = r'm:\deldir\1000x500000.bgen'
        #filepath = r'D:\OneDrive\Shares\bgenraaderpy\1x1000000.bgen'
        filepath = r'M:\del35\temp1024x16384-8.bgen'
        with open_bgen(filepath) as bgen2:
            print(bgen2.id[:5]) #other property arrays include risd,chrom,position,nallels, and allele_ids
            geno = bgen2.read(-1) # read the 200,000th variate's data
            #geno = bgen2.read() # read all, uses the ncombs from the first variant
            geno = bgen2.read(slice(5)) # read first 5, uses the ncombs from the first variant
            geno = bgen2.read(bgen2.chrom=='5',max_ncombs=4) # read chrom1, set max_combs explicitly
            geno, missing = bgen2.read(0,return_missing=True)
        bgen2 = open_bgen(filepath)
        geno, missing = bgen2.read(0,return_missing=True)
    if True:
        filepath = r'm:\deldir\2500x500000.bgen'
        with open_bgen(filepath,verbose=True) as bgen2:
            print(bgen2.read(0)[0,0,:])
            print(bgen2.read(-1)[0,0,:])
    print('!!!done')
