from pathlib import Path
from EnteroTyper.bin.enterobase_typer import type_sample


def get_sample_name_dict(indir: Path) -> dict:
    fasta_files = list(indir.glob("*.fna"))
    fasta_files += list(indir.glob("*.fasta"))
    fasta_files += list(indir.glob("*.fa"))

    sample_name_dict = {}
    for f in fasta_files:
        sample_name = f.with_suffix("").name
        sample_name_dict[sample_name] = f
    return sample_name_dict


def bulk_sample_typing(indir: Path, outdir: Path, database: Path) -> list:
    sample_name_dict = get_sample_name_dict(indir=indir)
    detailed_report_list = []
    for sample_name, assembly in sample_name_dict.items():
        sample_out_dir = outdir / sample_name
        detailed_report = type_sample(input_assembly=assembly, database=database, outdir=sample_out_dir,
                                      create_db=False, sample_name=sample_name)
        detailed_report_list.append(detailed_report)
    return detailed_report_list
