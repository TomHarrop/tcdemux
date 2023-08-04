#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import tempfile

def get_pool_barcode_seq(wildcards):
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

# config
sample_data_file = Path('data', 'samples.csv')
read_directory = Path('data', 'raw_data')
outdir = Path('output')
logdir = Path(outdir, 'logs')
keep_intermediate_files = True      # make this a cli option
workdir = outdir if keep_intermediate_files else tempfile.TemporaryDirectory()

# containers
bbmap = 'docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1'
cutadapt = 'docker://quay.io/biocontainers/cutadapt:4.4--py310h4b81fae_1'

# read the sample data and get a list of samples
sample_data = pd.read_csv(
    sample_data_file, index_col='name')

all_samples = sorted(set(sample_data.index))

# Check for the pool_name field.
# If the pool_name field is not present, we don't need to demultiplex internal
# barcodes, we only need to check the external barcodes for errors.
# Unfortunately that step is called demuxbyname.sh.
# Could do this with a conditional in the target, which would make the subsequent
# rules simpler.
samples_have_internal_barcodes = False
all_pools = []

if 'pool_name' in sample_data:
    samples_have_internal_barcodes = True
    all_pools = sorted(set(sample_data['pool_name']))

# all_samples = all_samples[0:1]

if samples_have_internal_barcodes:
    rule target:
        input:
            expand(
                Path(
                    workdir,
                    '020_demultiplex',
                    '{sample}_r1.fastq'
                    ),
                sample=all_samples
                )
else:
    rule target:
        input:
            expand(
                Path(
                    workdir,
                    '010_barcode-check',
                    '{sample}_r1.fastq'
                    ),
                sample=all_samples
                )

# demux using cutadapt
for mypool in all_pools:
    pool_sd = sample_data[sample_data['pool_name'] == mypool]
    pool_samples = sorted(set(pool_sd.index))
    rule:
        input:
            r1 = Path(
                workdir,
                '010_barcode-check',
                f'{mypool}_r1.fastq'
                ),
            r2 = Path(
                workdir,
                '010_barcode-check',
                f'{mypool}_r2.fastq'
                ),
            barcodes = Path(
                workdir,
                '020_demultiplex',
                f'{mypool}.barcodes.fasta'
                )
        output:
            r1 = expand(
                Path(
                    workdir,
                    '020_demultiplex',
                    '{sample}_r1.fastq'
                    ),
                sample = pool_samples
                ),
            r2 = expand(
                Path(
                    workdir,
                    '020_demultiplex',
                    '{sample}_r2.fastq'
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
            '-o {params.outdir}/{{name}}_r1.fastq '
            '-p {params.outdir}/{{name}}_r2.fastq '
            '{input.r1} '
            '{input.r2} '
            '&> {log}'
    rule:
        input:
            sample_data = sample_data_file
        output:
            barcode_file = Path(
                workdir,
                '020_demultiplex',
                f'{mypool}.barcodes.fasta'
                )
        params:
            pool = mypool
        script:
            'src/write_barcode_file.py'

# check barcodes in pooled fastq files
with tempfile.TemporaryDirectory() as tmpdir:
    rule check_pool_barcodes:
        input:
            unpack(get_pool_files)
        output:
            r1 = Path(
                workdir,
                '010_barcode-check',
                '{pool}_r1.fastq'
                ),
            r2 = Path(
                workdir,
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
                get_pool_barcode_seq(wildcards),
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
            '2> {log} '
            '&& '
            'mv {tmpdir}/{params.barcode_seq}_r1.fastq '
            '{output.r1} '
            '&& '
            'mv {tmpdir}/{params.barcode_seq}_r2.fastq '
            '{output.r2}'


# check barcodes in individual barcode files
def get_sample_barcode_seq(wildcards):
    sample = wildcards.sample
    my_sample_data = sample_data.loc[sample]
    i5_seq = my_sample_data['i5_index']
    i7_seq = my_sample_data['i7_index']
    i5_rc = Seq(i5_seq).reverse_complement()
    return(f'{i7_seq}+{i5_rc}')

def get_sample_files(wildcards):
    sample = wildcards.sample
    my_sample_data = sample_data.loc[sample]
    r1_file = my_sample_data['r1_file']
    r2_file = my_sample_data['r2_file']
    return {
        'r1': Path(read_directory, r1_file),
        'r2': Path(read_directory, r2_file)
    }

with tempfile.TemporaryDirectory() as tmpdir:
    rule check_sample_barcodes:
        input:
            unpack(get_sample_files)
        output:
            r1 = Path(
                workdir,
                '010_barcode-check',
                '{sample}_r1.fastq'
                ),
            r2 = Path(
                workdir,
                '010_barcode-check',
                '{sample}_r2.fastq'
                ),
            stats = Path(
                outdir,
                '010_barcode-check',
                '{sample}.demux_stats.txt'
                ),
        params:
            barcode_seq = lambda wildcards:
                get_sample_barcode_seq(wildcards),
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
                'check_pool_barcodes.{sample}.log'
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
            '2> {log} '
            '&& '
            'mv {tmpdir}/{params.barcode_seq}_r1.fastq '
            '{output.r1} '
            '&& '
            'mv {tmpdir}/{params.barcode_seq}_r2.fastq '
            '{output.r2}'