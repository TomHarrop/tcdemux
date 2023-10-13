#!/usr/bin/env python3

from pathlib import Path

version = config["version"]
testdir = Path(config["outdir"])

tcdemux = f"docker://quay.io/biocontainers/tcdemux:{version}"


rule target:
    input:
        Path(
            testdir,
            "pooled",
            "134473_LibID134573_GAP_BRF_H5TT7DRX3_ATCCACTG-ACGCACCT_S72_L002.r1.fastq.gz",
        ),
        Path(
            testdir,
            "unpooled",
            "134473_LibID134573_GAP_BRF_H5TT7DRX3_ATCCACTG-ACGCACCT_S72_L002.r1.fastq.gz",
        ),
        Path(testdir, "dirty.log"),
        Path(
            testdir,
            "missing",
            "this_sample_exists.r1.fastq.gz",
        ),


rule pooled:
    input:
        sd=Path("data/samples.csv"),
        adaptors=[
            Path("data/adaptors/alicia_adapters.fa"),
            Path("data/adaptors/TruSeq3-PE-2.fa"),
            Path("data/adaptors/bbmap_39.01_adaptors.fa"),
        ],
    output:
        f1=Path(
            testdir,
            "pooled",
            "134473_LibID134573_GAP_BRF_H5TT7DRX3_ATCCACTG-ACGCACCT_S72_L002.r1.fastq.gz",
        ),
    log:
        Path(testdir, "pooled.log"),
    params:
        readdir=Path("data/test_data"),
        outdir=lambda wildcards, output: Path(output.f1).parent,
        mem_gb=lambda wildcards, resources: int(resources.mem_mb // 1e3),
    threads: workflow.cores
    resources:
        mem_mb=int(16e3),
    container:
        tcdemux
    shell:
        "tcdemux "
        "--sample_data {input.sd} "
        "--adaptors {input.adaptors} "
        "--outdir {params.outdir} "
        "--read_directory {params.readdir} "
        "--threads {threads} "
        "--mem_gb {params.mem_gb} "
        "&> {log}"


rule unpooled:
    input:
        sd=Path("data/unpooled_samples.csv"),
        adaptors=[
            Path("data/adaptors/alicia_adapters.fa"),
            Path("data/adaptors/TruSeq3-PE-2.fa"),
            Path("data/adaptors/bbmap_39.01_adaptors.fa"),
        ],
    output:
        f1=Path(
            testdir,
            "unpooled",
            "134473_LibID134573_GAP_BRF_H5TT7DRX3_ATCCACTG-ACGCACCT_S72_L002.r1.fastq.gz",
        ),
    log:
        Path(testdir, "unpooled.log"),
    params:
        readdir=Path("data/test_data"),
        outdir=lambda wildcards, output: Path(output.f1).parent,
        mem_gb=lambda wildcards, resources: int(resources.mem_mb // 1e3),
    threads: workflow.cores
    resources:
        mem_mb=int(16e3),
    container:
        tcdemux
    shell:
        "tcdemux "
        "--sample_data {input.sd} "
        "--adaptors {input.adaptors} "
        "--outdir {params.outdir} "
        "--read_directory {params.readdir} "
        "--threads {threads} "
        "--mem_gb {params.mem_gb} "
        "&> {log}"


rule dirty:
    input:
        sd=Path("data/dirty_samples.csv"),
        adaptors=[
            Path("data/adaptors/alicia_adapters.fa"),
            Path("data/adaptors/TruSeq3-PE-2.fa"),
            Path("data/adaptors/bbmap_39.01_adaptors.fa"),
        ],
    log:
        Path(testdir, "dirty.log"),
    params:
        readdir=Path("data/test_data"),
        outdir=Path(testdir, "dirty"),
        mem_gb=lambda wildcards, resources: int(resources.mem_mb // 1e3),
    threads: workflow.cores
    resources:
        mem_mb=int(16e3),
    container:
        tcdemux
    shell:
        "tcdemux "
        "--sample_data {input.sd} "
        "--adaptors {input.adaptors} "
        "--outdir {params.outdir} "
        "--read_directory {params.readdir} "
        "--threads {threads} "
        "--mem_gb {params.mem_gb} "
        "&> {log}"


rule missing:
    input:
        sd=Path("data/samples_missing.csv"),
        adaptors=[
            Path("data/adaptors/alicia_adapters.fa"),
            Path("data/adaptors/TruSeq3-PE-2.fa"),
            Path("data/adaptors/bbmap_39.01_adaptors.fa"),
        ],
    output:
        f1=Path(
            testdir,
            "missing",
            "this_sample_exists.r1.fastq.gz",
        ),
    log:
        Path(testdir, "missing.log"),
    params:
        readdir=Path("data/test_data"),
        outdir=lambda wildcards, output: Path(output.f1).parent,
        mem_gb=lambda wildcards, resources: int(resources.mem_mb // 1e3),
    threads: workflow.cores
    resources:
        mem_mb=int(16e3),
    container:
        tcdemux
    shell:
        "tcdemux "
        "--sample_data {input.sd} "
        "--adaptors {input.adaptors} "
        "--outdir {params.outdir} "
        "--read_directory {params.readdir} "
        "--threads {threads} "
        "--mem_gb {params.mem_gb} "
        "&> {log}"
