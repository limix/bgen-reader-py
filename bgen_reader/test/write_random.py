import math
import os
import subprocess
from pathlib import Path

import numpy as np

from bgen_reader._helper import _log_in_place


def _write_random(
    filepath,
    nsamples,
    nvariants,
    bits=16,
    compression=None,
    chrom_count=22,
    verbose=True,
    seed=0,
    block_size=None,
    qctool_path=None,
    cleanup_temp_files=True,
):
    """
    If the external qctools is available, generates a random (unphased, diploid) BGEN file of any size.
    """
    qctool_path = qctool_path or os.environ.get("QCTOOLPATH")
    assert (
        qctool_path is not None
    ), "Bgen.write() requires a path to an external qctool program either via the qctool_path input or by setting the QCTOOLPATH environment variable."
    filepath = Path(filepath)
    gen_temp_filepath = filepath.with_suffix(".gen.temp")
    gen_filepath = filepath.with_suffix(".gen")
    sample_filepath = filepath.with_suffix(".sample")
    metadata2_filepath = filepath.with_suffix(".metadata.npz")

    # https://www.cog-genomics.org/plink2/formats#gen
    # https://web.archive.org/web/20181010160322/http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html
    # We need the +1 so that all three values will have enough precision to be very near 1
    # The max(3,..) is needed to even 1 bit will have enough precision in the gen file
    decimal_places = max(3, math.ceil(math.log(2 ** bits, 10)) + 1)

    def _format_function(num):
        return ("{0:." + str(decimal_places) + "f}").format(num)

    block_size = block_size or max((100 * 1000) // max(1, nsamples), 1)

    random_state = np.random.RandomState(seed)
    missing_rate = 0.218
    chrom_size = np.array(
        [
            263,
            255,
            214,
            203,
            194,
            183,
            171,
            155,
            145,
            144,
            144,
            143,
            114,
            109,
            106,
            98,
            92,
            85,
            67,
            72,
            50,
            56,
        ],
        dtype=np.int64,
    ) * int(
        1e6
    )  # The approximate size of human chromosomes in base pairs

    chromosomes = np.zeros((nvariants), dtype="int")
    positions = np.zeros((nvariants), dtype="int")

    chrom_total = chrom_size[:chrom_count].sum()
    chrom_step = chrom_total // nsamples
    chrom_start = 0
    chrom_size_so_far = 0
    for chrom_index in range(chrom_count):
        chrom_size_so_far += chrom_size[chrom_index]
        chrom_stop = chrom_size_so_far * nvariants // chrom_total
        chromosomes[chrom_start:chrom_stop] = chrom_index + 1
        random_increment = random_state.randint(chrom_step)
        positions[chrom_start:chrom_stop] = (
            np.arange(0, chrom_stop - chrom_start) * chrom_step + random_increment + 1
        )
        chrom_start = chrom_stop

    start = 0
    updater_freq = 10000
    index = -1
    with _log_in_place("_write_random", verbose) as updater:
        with gen_temp_filepath.open("w", newline="\n") as genfp:
            while start < nvariants:
                val = random_state.random(
                    (nsamples, min(block_size, nvariants - start), 3)
                )
                val /= val.sum(axis=2, keepdims=True)  # make probabilities sum to 1
                missing = random_state.rand(val.shape[0], val.shape[1]) < missing_rate
                val[missing, :] = np.nan

                for ivariants_in_block in range(val.shape[1]):
                    ivariants = start + ivariants_in_block
                    id = "SNP{0}".format(ivariants + 1)
                    rsid = "RS{0}".format(ivariants + 1)
                    chrom = chromosomes[ivariants]
                    pos = positions[ivariants]
                    genfp.write("{0} {1} {2} {3} A G".format(chrom, id, rsid, pos))
                    for isamples in range(nsamples):
                        index += 1
                        if updater_freq > 1 and index > 0 and index % updater_freq == 0:
                            updater(
                                "{0:,} of {1:,} ({2:2}%)".format(
                                    index + 1,
                                    nsamples * nvariants,
                                    100.0 * index / (nsamples * nvariants),
                                )
                            )
                        prob_dist = val[isamples, ivariants_in_block, :]
                        if not np.isnan(prob_dist).any():
                            s = " " + " ".join(
                                (_format_function(num) for num in prob_dist)
                            )
                            genfp.write(s)
                        else:
                            genfp.write(" 0 0 0")
                    genfp.write("\n")
                start += val.shape[1]
    # https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html
    with sample_filepath.open("w", newline="\n") as samplefp:
        samplefp.write("ID\n")
        samplefp.write("0\n")
        for isamples in range(nsamples):
            samplefp.write("sample_{0}\n".format(isamples))

    if gen_filepath.exists():
        gen_filepath.unlink()
    gen_temp_filepath.rename(gen_filepath)

    olddir = (
        os.getcwd()
    )  # This working directory stuff allows Windows to call qctools on a Linux subsystem
    os.chdir(filepath.parent)

    if filepath.exists():
        filepath.unlink()
    if metadata2_filepath.exists():
        metadata2_filepath.unlink()
    cmd = "{0} -g {1} -s {2} -og {3}{4}{5}".format(
        qctool_path,
        gen_filepath.name,
        sample_filepath.name,
        filepath.name,
        " -bgen-bits {0}".format(bits) if bits is not None else "",
        " -bgen-compression {0}".format(compression) if compression is not None else "",
    )
    if verbose:
        print("Calling qctool to convert *.gen to *.bgen")
    try:
        subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True, universal_newlines=True
        )
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)
        raise Exception("qctool command failed")
    if cleanup_temp_files:
        gen_filepath.unlink()
        sample_filepath.unlink()
    os.chdir(olddir)
