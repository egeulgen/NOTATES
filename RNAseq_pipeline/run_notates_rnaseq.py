#!/usr/bin/env python

##################################################
## Project: NOTATES RNAseq Pipeline
## Purpose: Wrapper script for the NOTATES RNAseq pipeline
## Date: Mar 2026
## Author: Ege Ulgen
##################################################

import subprocess
import argparse
import yaml
from typing import Any
from pathlib import Path

DEFAULT_STAR_ALIGN_OPTIONS = {
    "outSAMtype": "BAM Unsorted",
    "twopassMode": "Basic",
    "chimSegmentMin": 10,
    "chimJunctionOverhangMin": 10,
    "chimScoreMin": 1,
    "chimScoreDropMax": 30,
    "chimScoreJunctionNonGTAG": 0,
    "chimScoreSeparation": 1,
    "chimSegmentReadGapMax": 3,
    "alignSJstitchMismatchNmax": "5 -1 5 5",
    "chimOutType": "WithinBAM SoftClip",
}

DEFAULT_STAR_INDEX_OPTIONS = {
    "genomeSAsparseD": 3,
    "genomeSAindexNbases": 12,
}

CONFIG_FILE_NAME = "config.yaml"
DAG_FILE_NAME = "pipeline_dag.pdf"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="NOT - RNAseq-Pipeline",
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=40),
    )

    parser.add_argument(
        "-d",
        "--data-dir",
        type=Path,
        required=True,
        help="/path/to/raw/fastq/data",
    )

    parser.add_argument(
        "-i",
        "--sample-id",
        type=str,
        required=True,
        help="sample id string",
    )

    parser.add_argument(
        "-O",
        "--output_target",
        type=str,
        required=False,
        choices=["fastqc", "trimming", "quantification", "all"],
        default="all",
        help="Final output target; one of all (default), fastqc, trimming, mapping, quantification",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory",
    )

    parser.add_argument(
        "-R",
        "--reference-genome",
        type=Path,
        required=True,
        help="/path/to/reference/genome/fasta",
    )
    parser.add_argument(
        "-g",
        "--gtf-file",
        type=Path,
        required=True,
        help="/path/to/GTF/file",
    )
    parser.add_argument(
        "-b",
        "--fusion-blacklist",
        type=Path,
        required=True,
        help="/path/to/blacklist/file",
    )
    parser.add_argument(
        "-k",
        "--known-fusions",
        type=Path,
        required=True,
        help="/path/to/knowns/fusions/file",
    )

    parser.add_argument(
        "--star-index",
        type=Path,
        required=True,
        help="Path specifying where to create the STAR genome index (if it already exists, the pipeline automatically uses the index provided by the path)",
    )

    parser.add_argument(
        "--rsem-index",
        type=Path,
        required=True,
        help="Path specifying where to create the RSEM index (if it already exists, the pipeline automatically uses the index provided by the path)",
    )

    parser.add_argument(
        "-T",
        "--threads",
        type=int,
        required=False,
        default=8,
        help="Number of threads (default = 8)",
    )

    parser.add_argument(
        "--trim-galore-options",
        type=Path,
        required=False,
        help="Path to YAML config containing (optional) Trimgalore options",
    )

    parser.add_argument(
        "--star-align-options",
        type=Path,
        required=False,
        help="Path to YAML config containing (optional) STAR-align options",
    )

    parser.add_argument(
        "--star-index-options",
        type=Path,
        required=False,
        help="Path to YAML config containing (optional) STAR-index options",
    )

    parser.add_argument(
        "--create-DAG",
        action="store_true",
        help="specify to create the Directed Acyclic Graph of the workflow and exit (if not specified, pipeline runs)",
    )
    return parser.parse_args()


def load_options_from_yaml(opts_path: Path | None) -> dict[str, Any]:
    if opts_path and opts_path.exists():
        with opts_path.open("rt") as opt_handle:
            return yaml.safe_load(opt_handle)
    return {}


def generate_config_file(args: argparse.Namespace) -> Path:
    output_dir = args.output_dir
    config_content = {
        "sample_id": args.sample_id,
        "data_dir": str(args.data_dir.resolve()),
        "output_target": args.output_target,
        "output_dir": str(output_dir.resolve()),
        "reference_genome": str(args.reference_genome.resolve()),
        "reference_gtf": str(args.gtf_file.resolve()),
        "star_index": str(args.star_index.resolve()),
        "rsem_index": str(args.rsem_index.resolve()),
        "fusion_blacklist": str(args.fusion_blacklist.resolve()),
        "known_fusions": str(args.known_fusions.resolve()),
        "threads": args.threads,
        "trim_galore_options": load_options_from_yaml(args.trim_galore_options),
        "star_align_options": load_options_from_yaml(args.star_align_options)
        or DEFAULT_STAR_ALIGN_OPTIONS,
        "star_index_options": load_options_from_yaml(args.star_index_options)
        or DEFAULT_STAR_INDEX_OPTIONS,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir.joinpath(CONFIG_FILE_NAME)
    with config_path.open("wt") as outfile:
        yaml.dump(config_content, outfile, default_flow_style=False)
    return config_path


def create_dag_and_exit(config_path: Path) -> None:
    snakemake = subprocess.Popen(
        ["snakemake", "--configfile", str(config_path.resolve()), "--dag"],
        stdout=subprocess.PIPE,
    )
    dag_path = config_path.parent.joinpath(DAG_FILE_NAME)
    with dag_path.open("wb") as fp:
        subprocess.run(
            ["dot", "-Tpdf"],
            stdin=snakemake.stdout,
            stdout=fp,
            check=True,
        )


def execute_pipeline(config_path: Path, threads: int) -> None:
    subprocess.run(
        [
            "snakemake",
            "--configfile",
            str(config_path.resolve()),
            "--cores",
            str(threads),
        ],
        check=True,
    )


def main() -> None:
    arguments = parse_arguments()
    config_path = generate_config_file(arguments)

    if arguments.create_DAG:
        create_dag_and_exit(config_path)
    else:
        execute_pipeline(config_path, arguments.threads)


if __name__ == "__main__":
    main()
