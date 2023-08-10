#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import tempfile


def get_repair_input(wildcards):
    if samples_have_internal_barcodes:
        return {
        'r1': Path(
            workingdir,
            '020_demultiplex',
            '{sample}_r1.fastq'
            ),
        'r2': Path(
            workingdir,
            '020_demultiplex',
            '{sample}_r2.fastq'
            )
        }
    else:
        return{
        'r1': Path(
            workingdir,
            '010_barcode-check',
            '{sample}_r1.fastq'
            ),

        'r2': Path(
            workingdir,
            '010_barcode-check',
            '{sample}_r2.fastq'
            ),
        }

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
# read_directory = Path('data', 'test_data') # 1000 reads/file
outdir = Path('output')
logdir = Path(outdir, 'logs')
keep_intermediate_files = True      # make this a cli option
workingdir = outdir if keep_intermediate_files else tempfile.TemporaryDirectory()
adaptor_files = [
    'data/adaptors/alicia_adapters.fa',
    'data/adaptors/TruSeq3-PE-2.fa'
        ]
# adaptor references included with bbmap container
adaptor_path = '/usr/local/opt/bbmap-39.01-1/resources/adapters.fa'

# containers
bbmap = 'docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1'
cutadapt = 'docker://quay.io/biocontainers/cutadapt:4.4--py310h4b81fae_1'

# read the sample data and get a list of samples
sample_data = pd.read_csv(
    sample_data_file,
    index_col='name')

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

all_samples = all_samples[0:2]

rule target:
    input:
        expand(
            Path(
                outdir,
                '050_processed-files',
                '{sample}_r1.fastq.gz'
                ),
            sample=all_samples
            ),
        expand(
            Path(
                outdir,
                '050_processed-files',
                '{sample}.unpaired.fastq.gz'
                ),
            sample=all_samples
            )

# mask low complexity regions
rule mask:
    input:
        pipe = Path(
            workingdir,
            '040_trim',
            '{sample}.paired.fastq'
            )
    output:
        r1 = Path(
            outdir,
            '050_processed-files',
            '{sample}_r1.fastq.gz'
            ),
        r2 = Path(
            outdir,
            '050_processed-files',
            '{sample}_r2.fastq.gz'
            )
    log:
        Path(
            logdir,
            'mask',
            '{sample}.log'
            ),
    resources:
        time = lambda wildcards, attempt: 10 * attempt
    container:
        bbmap
    shell:
        'bbduk.sh '
        '-Xmx{resources.mem_mb}m '
        'threads={threads} '
        'in={input.pipe} '
        'int=t '
        'out={output.r1} '
        'out2={output.r2} '
        'zl=9 '
        'entropy=0.8 '
        'entropywindow=20 '
        'entropymask=t '
        '2> {log} '

rule mask_unpaired:
    input:
        Path(
            workingdir,
            '040_trim',
            '{sample}.unpaired.fastq'
            )
    output:
        Path(
            outdir,
            '050_processed-files',
            '{sample}.unpaired.fastq.gz'
            ),
    log:
        Path(
            logdir,
            'mask_unpaired',
            '{sample}.log'
            ),
    resources:
        time = lambda wildcards, attempt: 10 * attempt
    container:
        bbmap
    shell:
        'bbduk.sh '
        '-Xmx{resources.mem_mb}m '
        'threads={threads} '
        'in={input} '
        'int=f '
        'out={output} '
        'zl=9 '
        'entropy=0.8 '
        'entropywindow=20 '
        'entropymask=t '
        '2> {log} '

# trim adaptors
rule trim:
    input:
        pipe = Path(
            workingdir,
            '030_repair',
            '{sample}.fastq'
            ),
        adaptors = Path(
            outdir, 
            '040_trim',
            'adaptors.fasta'
            ),
    output:
        pipe = pipe(
            Path(
            workingdir,
            '040_trim',
            '{sample}.paired.fastq'
            )
            ),
        unpaired = temp(
            Path(
            workingdir,
            '040_trim',
            '{sample}.unpaired.fastq'
            )
            ),
        stats = Path(
            outdir,
            '040_trim',
            '{sample}.stats'
            )
    log:
        Path(
            logdir,
            'trim',
            '{sample}.log'
            )
    resources:
        time = lambda wildcards, attempt: 10 * attempt
    container:
        bbmap
    shell:
        'bbduk.sh '
        '-Xmx{resources.mem_mb}m '
        'threads={threads} '
        'in={input.pipe} '
        'int=t '
        'out=stdout.fastq '
        'outs={output.unpaired} '
        'ref={input.adaptors} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.stats} '
        '>> {output.pipe} '
        '2> {log} '

rule combine_adaptors:
    input:
        adaptor_files
    output:
        adaptors = Path(
            outdir, 
            '040_trim',
            'adaptors.fasta'
            ),
        duplicates = Path(
            outdir, 
            '040_trim',
            'duplicated_adaptors.fasta'
            )
    log:
        Path(
            logdir,
            'combine_adaptors.log')
    threads:
        2
    resources:
        time = lambda wildcards, attempt: 10 * attempt
    container:
        bbmap
    shell:
        'cat {input} {adaptor_path} | '
        'dedupe.sh '
        '-Xmx{resources.mem_mb}m '
        'in=stdin.fasta '
        'out={output.adaptors} '
        'outd={output.duplicates} '
        'absorbcontainment=f '
        'absorbrc=f '
        'ascending=t '
        'exact=t '
        'maxedits=0 '
        'maxsubs=0 '
        'sort=name '
        'touppercase=t '
        'uniquenames=t '
        '2> {log} '


# double check pairing
rule repair:
    input:
        unpack(get_repair_input)
    output:
        pipe(
            Path(
            workingdir,
            '030_repair',
            '{sample}.fastq'
            )
            )
    log:
        Path(
            logdir,
            'repair',
            '{sample}.log'),
    resources:
        time = lambda wildcards, attempt: 10 * attempt
    container:
        bbmap
    shell:
        'repair.sh '
        '-Xmx{resources.mem_mb}m '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'outs=/dev/null '
        'repair=t '
        'tossbrokenreads=t '
        'tossjunk=t '
        '>> {output} '
        '2> {log}'

# demux using cutadapt
for mypool in all_pools:
    pool_sd = sample_data[sample_data['pool_name'] == mypool]
    pool_samples = sorted(set(pool_sd.index))
    rule:
        input:
            r1 = Path(
                workingdir,
                '010_barcode-check',
                f'{mypool}_r1.fastq'
                ),
            r2 = Path(
                workingdir,
                '010_barcode-check',
                f'{mypool}_r2.fastq'
                ),
            barcodes = Path(
                workingdir,
                '020_demultiplex',
                f'{mypool}.barcodes.fasta'
                )
        output:
            r1 = expand(
                Path(
                    workingdir,
                    '020_demultiplex',
                    '{sample}_r1.fastq'
                    ),
                sample = pool_samples
                ),
            r2 = expand(
                Path(
                    workingdir,
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
                'cutadapt',
                f'{mypool}.log'
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
                workingdir,
                '020_demultiplex',
                f'{mypool}.barcodes.fasta'
                )
        params:
            pool = mypool
        script:
            'src/write_barcode_file.py'

# check barcodes in pooled fastq files
with tempfile.TemporaryDirectory() as rule_tmpdir:
    rule check_pool_barcodes:
        input:
            unpack(get_pool_files)
        output:
            r1 = Path(
                workingdir,
                '010_barcode-check',
                '{pool}_r1.fastq'
                ),
            r2 = Path(
                workingdir,
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
                rule_tmpdir,
                '%_r1.fastq'
                ).as_posix(),
            out2 = Path(
                rule_tmpdir,
                '%_r2.fastq'
                ).as_posix()
        log:
            Path(
                logdir,
                'check_pool_barcodes',
                '{pool}.log'
                )
        threads:
            5
        resources:
            mem_gb = 24
        container:
            bbmap
        shell:
            'streams=$(( {threads}/2 )) ; '
            'mkdir -p {rule_tmpdir} ; '
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
            'mv {rule_tmpdir}/{params.barcode_seq}_r1.fastq '
            '{output.r1} '
            '&& '
            'mv {rule_tmpdir}/{params.barcode_seq}_r2.fastq '
            '{output.r2}'


# check barcodes in individual barcode files
with tempfile.TemporaryDirectory() as tmpdir:
    rule check_sample_barcodes:
        input:
            unpack(get_sample_files)
        output:
            r1 = Path(
                workingdir,
                '010_barcode-check',
                '{sample}_r1.fastq'
                ),
            r2 = Path(
                workingdir,
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
                'check_pool_barcodes',
                '{sample}.log'
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