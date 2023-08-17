#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools import find_packages


# load README.rst
def readme():
    with open('README.rst') as file:
        return file.read()

setup(
    name='tcdemux',
    version='0.0.7',
    description=(
        'python3 wrapper for demultiplexing '
        'target capture sequencing results'
        ),
    long_description=readme(),
    url='https://github.com/tomharrop/tcdemux',
    author='Tom Harrop',
    author_email='twharrop@gmail.com',
    license='GPL-3',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.81',
        'cutadapt>=4.4',
        'pandas>=2.0.3',
        'snakemake>=7.31.0'
    ],
    entry_points={
        'console_scripts': [
            'tcdemux = tcdemux.__main__:main'
            ],
    },
    scripts={
        'tcdemux/src/write_barcode_file.py',
    },
    package_data={
        'tcdemux': [
            'Snakefile',
            'README.rst'
        ],
    },
    zip_safe=False)
