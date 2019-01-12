from bgen_reader import example_files, read_bgen


def test_large_file(capsys):
    with capsys.disabled():
        with example_files("large.bgen") as filepath:
            data = read_bgen(filepath, verbose=True)
            variants = data["variants"]
            samples = data["samples"]
            genotype = data["genotype"]
            pass
