import click
import logging
from bin.enterobase_typer import type_sample, makeblastdb_database
from bin.sequence_concatenation_pipeline import sequence_concatenation_pipeline
from pathlib import Path


def get_sample_name_dict(indir: Path) -> dict:
    fasta_files = list(indir.glob("*.fna"))
    fasta_files += list(indir.glob("*.fasta"))
    fasta_files += list(indir.glob("*.fa"))

    sample_name_dict = {}
    for f in fasta_files:
        sample_name = f.with_suffix("").name
        sample_name_dict[sample_name] = f
    return sample_name_dict


if __name__ == "__main__":
    main()
