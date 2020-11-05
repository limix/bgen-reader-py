from ._environment import BGEN_READER_CACHE_HOME
from ._file import file_hash

_filenames = {
    "complex.23bits.bgen": "af9a0232fdf339512ef3d23c39a4c40dc8908e79c0aafae62063e6f3bb0c6fe3",
    "complex.23bits.no.samples.bgen": "25d30a4e489da1aeb05f9893af98e8bf3b09d74db2982bf1828f8c8565886fc6",
    "complex.bgen": "af9a0232fdf339512ef3d23c39a4c40dc8908e79c0aafae62063e6f3bb0c6fe3",
    "complex.sample": "19a9149e0551f2862c26be48e006b8ac8cd0bd9ca2793ca82ca4b63a1c16083f",
    "example.32bits.bgen": "76386a9f100a3f1bd5e88bb8ba6a81993c5673c0b39cef9fd786298adac2fbd5",
    "example.bgen": "76386a9f100a3f1bd5e88bb8ba6a81993c5673c0b39cef9fd786298adac2fbd5",
    "haplotypes.bgen": "84e0b59efcc83c7c305cf5446e5dc26b49b15aeb4157a9eb36451376ce3efe4c",
    "haplotypes.bgen.metadata.corrupted": "8f55628770c1ae8155c1ced2463f15df80d32bc272a470bb1d6b68225e1604b1",
    "haplotypes.bgen.metadata.valid": "68215ce14966a07742979f666b9dfea8f1becafc9a15185eb01009654c7f5fe0",
    "wrong.metadata": "f746345605150076f3234fbeea7c52e86bf95c9329b2f08e1e3e92a7918b98fb",
}


def example_filepath(filename: str):
    import requests

    url = "https://bgen-examples.s3.amazonaws.com"

    if filename not in _filenames:
        raise ValueError(f"Unknown filename {filename}.")

    test_data_folder = BGEN_READER_CACHE_HOME / "test_data"
    test_data_folder.mkdir(parents=True, exist_ok=True)
    filepath = test_data_folder / filename

    if filepath.exists() and file_hash(filepath) != _filenames[filename]:
        filepath.unlink()

    if not filepath.exists():
        r = requests.get(f"{url}/{filename}")
        r.raise_for_status()
        with open(filepath, "wb") as f:
            f.write(r.content)

    if file_hash(filepath) != _filenames[filename]:
        msg = (
            f"Hash mismatch:\n"
            f"  ACTUAL : {file_hash(filepath)}\n"
            f"  DESIRED: {_filenames[filename]}"
        )
        raise RuntimeError(msg)

    return filepath
