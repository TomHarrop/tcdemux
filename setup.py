#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools import find_packages


# load README.rst
def readme():
    with open("tcdemux/README.rst") as file:
        return file.read()


setup(
    name="tcdemux",
    version="0.1.1",
    description=(
        "python3 wrapper for demultiplexing target capture sequencing results"
    ),
    long_description=readme(),
    url="https://github.com/tomharrop/tcdemux",
    author="Tom Harrop",
    author_email="twharrop@gmail.com",
    license="GPL-3",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.81",
        "cutadapt>=4.5",
        "pandas>=2.1.1",
        "snakemake>=7.32.4,<8.0.0",
    ],
    entry_points={
        "console_scripts": ["tcdemux = tcdemux.__main__:main"],
    },
    scripts={
        "tcdemux/src/parse_cutadapt_stats.py",
        "tcdemux/src/plot_adaptor_stats.R",
        "tcdemux/src/plot_duplication_stats.R",
        "tcdemux/src/process_step_logs.py",
        "tcdemux/src/write_barcode_file.py",
    },
    package_data={
        "tcdemux": ["Snakefile", "README.rst"],
    },
    zip_safe=False,
)
