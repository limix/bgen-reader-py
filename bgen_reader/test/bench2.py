import numpy as np
import time
from pathlib import Path
from bgen_reader import open_bgen

def read_nth(nth, bgen2):

    bgen2_t0 = time.time()
    proball_bgen2 = bgen2.read(nth)
    #for ivariant in range(0,bgen2.nvariants)nth:
    #    prob1_bgen2 = bgen2.read(ivariant)#nth
    diff_bgen2 = time.time()-bgen2_t0
    print(f"bgen2: reading every {nth}th variant takes {diff_bgen2} seconds")
    prob1_bgen2 = proball_bgen2[:,-1,:]

    #cb_t0 = time.time()
    #proball_cb = np.empty((bgen2.nsamples,len(offset_list[nth]),3),order='F') # Much faster than order='C'
    ##    proball_cb = []
    #for ivariant, offset in enumerate(offset_list[nth]):
    #    prob1 = cb.read_probability(offset)
    #    proball_cb[:,ivariant,:] = prob1
    #    #proball_cb.append(prob1)
    ##proball_cb = np.array(proball_cb)
    #diff_cb = time.time()-cb_t0
    #print(f"cbgen: reading every {nth}th variant takes {diff_cb} seconds")
    #prob1_cb = proball_cb[:,-1,:]

    #assert np.allclose(prob1_bgen2.reshape(-1,prob1_bgen2.shape[-1]),prob1_cb,equal_nan=True)


#!!!cmk delete this whole file?
if __name__ == "__main__":
    #filename = "merged_487400x220000.bgen"
    filename = "1000x500000.bgen"

    cache_home = Path(r'M:\deldir\genbgen\good')
    filepath = cache_home / "test_data" / filename
    assert filepath.exists()
    bgen2 = open_bgen(filepath,verbose=True)

    #read_nth(np.s_[::5000], cb, offset_list, bgen2)
    #read_nth(np.s_[::500], cb, offset_list, bgen2)
    #read_nth(np.s_[::50], cb, offset_list, bgen2)
    read_nth(np.s_[::10], bgen2)
    #read_nth(np.s_[::5], cb, offset_list, bgen2)
    #read_nth(np.s_[::500], cb, offset_list, bgen2)
    #read_nth(np.s_[200*1000:200*1000+50], cb, offset_list, bgen2)
    #read_nth(np.s_[200*1000:200*1000+500], cb, offset_list, bgen2)



    del bgen2

    print("!!!cmk")