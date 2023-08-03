#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import tempfile

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
cutadapt = 'docker://quay.io/biocontainers/cutadapt:4.4--py310h4b81fae_1'

sample_data = pd.read_csv(
    sample_data_file, index_col='Name')

all_samples = sorted(set(sample_data.index))
all_pools = sorted(set(sample_data['pool_name']))

all_samples = all_samples[0:1]


rule target:
    input:
        expand(
            Path(
                outdir,
                '020_demultiplex',
                '{sample}_r1.fastq.gz'
                ),
            sample=all_samples
            )

for mypool in all_pools:
    pool_sd = sample_data[sample_data['pool_name'] == mypool]
    pool_samples = sorted(set(pool_sd.index))
    rule:
        input:
            r1 = Path(
                outdir,
                '010_barcode-check',
                f'{mypool}_r1.fastq'
                ),
            r2 = Path(
                outdir,
                '010_barcode-check',
                f'{mypool}_r2.fastq'
                ),
            barcodes = Path(
                outdir,
                '020_demultiplex',
                f'{mypool}.barcodes.fasta'
                )
        output:
            r1 = expand(
                Path(
                    outdir,
                    '020_demultiplex',
                    '{sample}_r1.fastq.gz'
                    ),
                sample = pool_samples
                ),
            r2 = expand(
                Path(
                    outdir,
                    '020_demultiplex',
                    '{sample}_r2.fastq.gz'
                    ),
                sample = pool_samples
                )
        params:
            outdir = lambda wildcards, output:
                Path(output[0]).parent.as_posix()
        log:
            Path(
                logdir,
                f'cutadapt.{mypool}.log'
                )
        threads:
            16
        container:
            cutadapt
        shell:
            'cutadapt '
            '-j {threads} '
            '-e 0 '
            '--no-indels '
            '--pair-adapters '
            '-g ^file:{input.barcodes} '
            '-G ^file:{input.barcodes} '
            '-o {params.outdir}/{{name}}_r1.fastq.gz '
            '-p {params.outdir}/{{name}}_r2.fastq.gz '
            '{input.r1} '
            '{input.r2} '
            '&> {log}'
    rule:
        input:
            sample_data = sample_data_file
        output:
            barcode_file = Path(
                outdir,
                '020_demultiplex',
                f'{mypool}.barcodes.fasta'
                )
        params:
            pool = mypool
        script:
            'src/write_barcode_file.py'


with tempfile.TemporaryDirectory() as tmpdir:
    rule check_pool_barcodes:
        input:
            unpack(get_pool_files)
        output:
            r1 = Path(
                outdir,
                '010_barcode-check',
                '{pool}_r1.fastq'
                ),
            r2 = Path(
                outdir,
                '010_barcode-check',
                '{pool}_r2.fastq'
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
                tmpdir,
                '%_r1.fastq'
                ).as_posix(),
            out2 = Path(
                tmpdir,
                '%_r2.fastq'
                ).as_posix()
        log:
            Path(
                logdir,
                'check_pool_barcodes.{pool}.log'
                )
        threads:
            5
        resources:
            mem_gb = 24
        container:
            bbmap
        shell:
            'streams=$(( {threads}/2 )) ; '
            'demuxbyname.sh ' 
            'delimiter=: prefixmode=f ' # use the last column
            'names={params.barcode_seq} '
            'out={params.out} '
            'out2={params.out2} '
            'outu=/dev/null '
            'stats={output.stats} '
            'in={input.r1} '
            'in2={input.r2} '
            'streams=$streams '
            '-Xmx{resources.mem_gb}g '
            # 'zl=9 '
            '2> {log} '
            '&& '
            'mv {tmpdir}/{params.barcode_seq}_r1.fastq '
            '{output.r1} '
            '&& '
            'mv {tmpdir}/{params.barcode_seq}_r2.fastq '
            '{output.r2}'
