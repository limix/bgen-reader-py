from __future__ import unicode_literals

import os
import sys

from setuptools import find_packages, setup

try:
    import pypandoc
    long_description = pypandoc.convert_file('README.md', 'rst')
except (OSError, IOError, ImportError):
    long_description = open('README.md').read()


def setup_package():
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
    pytest_runner = ['pytest-runner>=2.9'] if needs_pytest else []

    setup_requires = ['cffi>=1.8'] + pytest_runner
    install_requires = [
        'cffi>=1.8', 'tqdm>=4.14', 'dask[array,bag,dataframe,delayed]>=0.15',
        'scipy>=0.18', 'pandas>=0.19.2'
    ]
    tests_require = ['pytest', 'pytest-pep8']

    metadata = dict(
        name='bgen-reader',
        version='0.1.9',
        maintainer="Danilo Horta",
        maintainer_email="horta@ebi.ac.uk",
        author="Danilo Horta",
        author_email="horta@ebi.ac.uk",
        license="MIT",
        description="Bgen file format reader.",
        long_description=long_description,
        url='https://github.com/limix/bgen-reader-py',
        packages=find_packages(),
        zip_safe=False,
        install_requires=install_requires,
        setup_requires=setup_requires,
        tests_require=tests_require,
        include_package_data=True,
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        cffi_modules=["bgen_reader/_build.py:ffibuilder"])

    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)


if __name__ == '__main__':
    setup_package()
