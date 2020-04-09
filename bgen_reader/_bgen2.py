import os
import numpy as np
from contextlib import contextmanager
from tempfile import mkdtemp
import shutil
from pathlib import Path
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

#!!!cmk should this be 'read_bgen2' instead?
@contextmanager
def bgen_reader2(filename, sample=None, verbose=False):
    bgen2 = Bgen2(filename,sample=sample,verbose=verbose)
    yield bgen2
    del bgen2

#!!!cmk can we have it both a contextmanger and not?
#!!!cmk add type info
#!!!cmk write doc
#!!!cmk test doc
#!!!cmk create test cases that generate big datasets (if qctool is available)
#!!!cmk ok to have metadata2.npz location be fixed for now?
class Bgen2(object):
    def __init__(self, filename, sample=None, verbose=False):
        self._verbose = verbose
        self.filename = filename

        assert os.path.exists(filename), "Expect file to exist ('{0}')".format(filename) #!!!cmk still useful?

        self._bgen_context_manager = bgen_file(Path(filename))
        self._bgen = self._bgen_context_manager.__enter__()


        self.samples = np.array(self._get_samples(sample),dtype='str')

        metadata2 = self.filename + ".metadata2.npz"
        if os.path.exists(metadata2):
            d = np.load(metadata2) #!!!cmk could do memory mapping instead
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

    def _get_samples(self,sample_file):
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

    #!!!cmk add an nvariants property and nsamples
    #!!!cmk should have dtype (because float32 is often enough and is 1/2 the size) and order
    def read(self, variants=None, max_ncombs=None): #!!!cmk also allow samples to be selected?
        #!!!cmk allow single ints, lists of ints, lists of bools, None, and slices
        #!!!cmk could allow strings (variant names) and lists of strings

        max_ncombs = max_ncombs or self.max_ncombs

        if type(variants) is np.int: #!!!make sure this works with all int types
            variants = [variants]
        if variants is None:
            vaddr = self.vaddr
            ncombs = self.ncombs
        else:
            vaddr = self.vaddr[variants]
            ncombs = self.ncombs[variants]

        #allocating p only once make reading 10x5M data 30% faster
        val = np.full((len(self.samples), len(vaddr), max_ncombs), np.nan, order='F', dtype='float64') #!!!cmk test on selecting zero variants
        p = None

        #LATER multithread?
        #!!!cmk if verbose is true, give some status
        for out_index,vaddr0 in enumerate(vaddr):
            if p is None or ncombs[out_index] != p.shape[-1]:
                p = np.full((len(self.samples), ncombs[out_index]), np.nan, order='C', dtype='float64')
            genotype = lib.bgen_file_open_genotype(self._bgen._bgen_file, vaddr0)
            lib.bgen_genotype_read(genotype, ffi.cast("double *", p.ctypes.data))
            #ploidy = asarray([lib.bgen_ploidy(vg, i) for i in range(nsamples)], int) #!!!cmk what is this? It will likely be a different call
            #missing = asarray([lib.bgen_missing(vg, i) for i in range(nsamples)], bool) #!!!cmk why is this need instead of just seeing nan,nan,nan
            lib.bgen_genotype_close(genotype)
            val[:,out_index,:ncombs[out_index]] = p
        return val

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

    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self.filename)

    def __del__(self):
        if hasattr(self,'_bgen_context_manager') and self._bgen_context_manager is not None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self._bgen_context_manager.__exit__(None,None,None)

if __name__ == "__main__":
    if True:
        #filename = r'm:\deldir\1000x500000.bgen'
        #filename = r'D:\OneDrive\Shares\bgenreaderpy\1x1000000.bgen'
        filename = r'M:\del35\temp1024x16384-8.bgen'
        with bgen_reader2(filename) as bgen2:
            print(bgen2.id[:5]) #other property arrays include risd,chrom,position,nallels, and allele_ids
            geno = bgen2.read(-1) # read the 200,000th variate's data
            #geno = bgen2.read() # read all, uses the ncombs from the first variant
            geno = bgen2.read(slice(5)) # read first 5, uses the ncombs from the first variant
            #!!!cmk is there any prettier way for users to specify slices?
            geno = bgen2.read(bgen2.chrom=='5',max_ncombs=4) # read chrom1, set max_combs explicitly
    if True:
        filename = r'm:\deldir\2500x500000.bgen'
        with bgen_reader2(filename,verbose=True) as bgen2:
            print(bgen2.read(0)[0,0,:])
            print(bgen2.read(-1)[0,0,:])
    print('!!!done')
