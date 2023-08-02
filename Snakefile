#!/usr/bin/env python3

from Bio.Seq import Seq
from pathlib import Path
import pandas as pd

def get_barcode_seq(wildcards):
    pool = wildcards.pool
    pool_data = sample_data[
        sample_data['pool_name'] == pool
        ]
    i5_seq = pool_data['i5_index'][0]
    i7_seq = pool_data['i7_index'][0]
    i5_rc = Seq(i5_seq).reverse_complement()
    return(f'{i7_seq}+{i5_rc}')

def get_pool_files(wildcards):
    pool = wildcards.pool
    pool_data = sample_data[
        sample_data['pool_name'] == pool
        ]
    r1_file = pool_data['r1_file'][0]
    r2_file = pool_data['r2_file'][0]
    return {
        'r1': Path(read_directory, r1_file),
        'r2': Path(read_directory, r2_file)
    }


sample_data_file = Path('data', 'samples.csv')
read_directory = Path('data', 'raw_data')
outdir = Path('output')
logdir = Path(outdir, 'logs')

# containers
bbmap = 'docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1'

sample_data = pd.read_csv(
    sample_data_file, index_col='Name')

all_samples = sorted(set(sample_data.index))
all_pools = sorted(set(sample_data['pool_name']))


rule target:
    input:
        expand(
            Path(
                outdir,
                '010_barcode-check',
                '{pool}.demux_stats.txt'
                ),
            pool = all_pools
            )

rule check_pool_barcodes:
    input:
        unpack(get_pool_files)
    output:
        r1 = Path(
            outdir,
            '010_barcode-check',
            '{pool}_r1.fastq.gz'
            ),
        r2 = Path(
            outdir,
            '010_barcode-check',
            '{pool}_r2.fastq.gz'
            ),
        stats = Path(
            outdir,
            '010_barcode-check',
            '{pool}.demux_stats.txt'
            ),
    params:
        barcode_seq = lambda wildcards:
            get_barcode_seq(wildcards),
        out = Path(
            outdir,
            'tmp/reads/%_r1.fastq'
            ).as_posix(),
        out2 = Path(
            outdir,
            'tmp/reads/%_r2.fastq'
            ).as_posix()
    log:
        Path(
            logdir,
            'check_pool_barcodes.{pool}.log'
            )
    container:
        bbmap
    shell:
        'demuxbyname.sh ' 
        'delimiter=: prefixmode=f ' # use the last column
        'names={params.barcode_seq} '
        'out={params.out} '
        'out2={params.out2} '
        'in={input.r1} '
        'in2={input.r2} '
        '2> {log}'




