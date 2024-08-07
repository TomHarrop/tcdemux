#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from snakemake.logging import logger
import argparse
from importlib.metadata import files, version
import psutil
import snakemake

#############
# FUNCTIONS #
#############


def parse_arguments():
    parser = argparse.ArgumentParser(prog="tcdemux")
    parser.add_argument("-n", help="Dry run", dest="dry_run", action="store_true")
    default_threads = 4
    parser.add_argument(
        "--threads",
        help=("Number of threads. Default: %i" % default_threads),
        metavar="int",
        type=int,
        dest="threads",
        default=default_threads,
    )
    default_mem_gb = int(psutil.virtual_memory().available * 0.8 // 1e9)
    parser.add_argument(
        "--mem_gb",
        help=("Amount of RAM in GB. Default: %i" % default_mem_gb),
        metavar="int",
        type=int,
        dest="mem_gb",
        default=default_mem_gb,
    )
    parser.add_argument(
        "--restart_times",
        required=False,
        help="number of times to restart failing jobs (default 0)",
        type=int,
        dest="restart_times",
        default=0,
    )
    parser.add_argument(
        "--sample_data",
        required=True,
        help="Sample csv (see README)",
        type=str,
        dest="sample_data_file",
    )
    parser.add_argument(
        "--read_directory",
        required=True,
        help="Directory containing the read files",
        type=str,
        dest="read_directory",
    )
    parser.add_argument(
        "--adaptors",
        required=True,
        help="FASTA file(s) of adaptors. Multiple adaptor files can be used.",
        type=str,
        dest="adaptor_files",
        nargs="+",
    )
    parser.add_argument(
        "--outdir", required=True, help="Output directory", type=str, dest="outdir"
    )
    parser.add_argument(
        "--keep_intermediate_files", action=argparse.BooleanOptionalAction
    )
    parser.add_argument(
        "--qtrim",
        help="Trim right end of reads to remove bases with quality below trimq.",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--trimq",
        type=float,
        dest="trimq",
        help="Regions with average quality BELOW this will be trimmed, if qtrim is enabled",
        default=6.0,
    )

    args = vars(parser.parse_args())
    return args


########
# MAIN #
########


def main():
    # get the tcdemux version
    pkg_version = version(__package__)
    logger.info(f"tcdemux version {pkg_version}")

    # get the snakefile
    snakefile = [p.locate() for p in files(__package__) if "Snakefile" in str(p)][0]
    logger.debug(f"Using snakefile {snakefile}")

    # get args
    args = parse_arguments()
    logger.debug(f"Entrypoint args\n{args}")

    # define resources
    smk_resources = {"mem_mb": int(args["mem_gb"] * 1e3)}

    # run the pipeline
    snakemake.snakemake(
        snakefile=snakefile,
        config=args,
        cores=args["threads"],
        resources=smk_resources,
        overwrite_resource_scopes={
            "mem_gb": "global",
            "threads": "global",
        },
        printshellcmds=logger.printshellcmds,
        dryrun=True if args["dry_run"] else False,
        restart_times=args["restart_times"],
        lock=False,
    )
