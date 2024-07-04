#!/usr/bin/env python3

from pathlib import Path

version = config["version"]
testdir = Path(config["outdir"])
datadir = Path("test-data")
adaptors = (
    [
        Path(datadir, "adaptors/our_adaptors.fa"),
        Path(datadir, "adaptors/TruSeq3-PE-2.fa"),
        Path(datadir, "adaptors/bbmap_39.01_adaptors.fa"),
    ],
)


tcdemux = f"docker://quay.io/biocontainers/tcdemux:{version}"


rule target:
    input:
        expand(
            Path(
                testdir,
                "{testname}",
                "{trim}",
                "done",
            ),
            testname=["pooled", "unpooled", "missing"],
            trim=["qtrim", "no-qtrim"],
        ),
        Path(testdir, "dirty.log"),


rule tcdemux:
    input:
        sd=Path(datadir, "samples_{testname}.csv"),
        adaptors=adaptors,
    output:
        flag=touch(
            Path(
                testdir,
                "{testname}",
                "{trim}",
                "done",
            )
        ),
    log:
        Path(testdir, "{testname}.{trim}.log"),
    params:
        readdir=Path(datadir, "raw_reads"),
        outdir=lambda wildcards, output: Path(output.flag).parent,
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
        "--{wildcards.trim} "
        "&> {log}"


rule dirty:
    input:
        sd=Path(datadir, "samples_dirty.csv"),
        adaptors=adaptors,
    log:
        Path(testdir, "dirty.log"),
    params:
        readdir=Path(datadir, "raw_reads"),
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
