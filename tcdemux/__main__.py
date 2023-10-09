#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pkg_resources import resource_filename
from snakemake.logging import logger
import argparse
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

    args = vars(parser.parse_args())
    return args


########
# MAIN #
########


def main():
    # get the snakefile
    snakefile = resource_filename(__name__, "Snakefile")
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
        overwrite_resource_scopes={"mem_gb": "global"},
        printshellcmds=logger.printshellcmds,
        dryrun=True if args["dry_run"] else False,
        restart_times=args["restart_times"],
        lock=False,
    )
