import os
import click
import logging
import pandas as pd
import multiprocessing
from pathlib import Path
from subprocess import Popen
from enterobase_typer import get_database_files, run_subprocess, create_outdir, generate_cgmlst_report
from sequence_type_comparison import call_sequence_comparison


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    value = Path(value)
    return value


@click.command(help="Takes a list of target *.BLASTn_Detailed_Report.tsv files, followed by several options. "
                    "Extracts sequences from the BLASTn report files as .fasta files, aligns them all with MUSCLE, and "
                    "finally runs FASconCAT. ")
@click.option('-o', '--out_dir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Root directory to store all output files',
              callback=convert_to_path)
@click.option('-db', '--database',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Path to your MLST database',
              callback=convert_to_path)
@click.option('-v', '--verbose',
              is_flag=True,
              default=False,  # Set this to false eventually
              help='Set this flag to enable more verbose logging.')
@click.argument('targets', nargs=-1, type=click.Path(exists=True))
def main(targets, out_dir, database, verbose):
    if verbose:
        logging.basicConfig(
            format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
            level=logging.DEBUG,
            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(
            format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')

    fasconcat_exec = Path(__file__).parent / 'FASconCAT-G_v1.04.pl'
    if not fasconcat_exec.exists():
        logging.error(f"ERROR: Could not find {fasconcat_exec}")
        quit()

    logging.info("Started FASconCAT pipeline")

    fasconcat_pipeline(targets=targets, database=database, out_dir=out_dir, fasconcat_exec=fasconcat_exec)

    logging.info("Script complete")


def fasconcat_pipeline(targets: list, database: Path, out_dir: Path, fasconcat_exec: Path):
    targets = [Path(target) for target in targets]
    database_files = get_database_files(database=database)
    create_outdir(out_dir=out_dir)

    report_dict = {}
    for i, target in enumerate(targets):
        target_name = target.with_suffix("").name + f'_{i}'
        report_dict[target_name] = target

    logging.info("Generating .FASTA files")
    df_params = [(database_file, out_dir, report_dict) for database_file in database_files]
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as p:
        p.starmap(write_fasta, df_params)

    logging.info("Aligning .FASTA files with MUSCLE")
    fasta_files = list(out_dir.glob("*.fasta"))
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as p:
        p.map(call_muscle, fasta_files)

    # Cleanup
    logging.info("Deleting interim files")
    [os.remove(str(fasta)) for fasta in fasta_files]

    # Grab newly aligned fasta files
    aligned_fasta_files = list(out_dir.glob("*.align.fasta"))

    # Create multifasta with fasconcat
    call_fasconcat(target_dir=out_dir, fasconcat_exec=fasconcat_exec)

    # Remove remaining fasta files
    [os.remove(str(fasta)) for fasta in aligned_fasta_files]

    # Sequence type comparison
    sequence_type_report_list = []
    with click.progressbar(targets, length=len(targets), label='Extracting sequence types') as bar:
        for target in bar:
            df = pd.read_csv(target, sep='\t')
            cgmlst_allele_report = generate_cgmlst_report(df=df, out_dir=out_dir,
                                                          sample_name=target.name.split("_")[0])
            sequence_type_report_list.append(cgmlst_allele_report)

    call_sequence_comparison(targets=sequence_type_report_list, out_dir=out_dir / 'sequence_type_comparisons')

    # Cleanup reports
    [os.remove(str(report)) for report in sequence_type_report_list]


def get_top_qseq(df: pd.DataFrame, locus: str) -> str:
    locus_df = df[df['locus'] == locus]
    qseq = locus_df['qseq_strand_aware'].iloc[0]
    if qseq == "NA":
        return ""
    return qseq


def read_blast_report(report: Path):
    df = pd.read_csv(report, delimiter="\t")
    return df


def write_fasta(database_file: Path, outdir: Path, report_dict: dict):
    locus = database_file.with_suffix("").name
    locus_length = 0
    logging.debug(f"{locus}...")
    outfile = str(outdir / database_file.with_suffix(".fasta").name)
    with open(outfile, "w") as out:
        for target, report_path in report_dict.items():
            out.write(f">{target}\n")
            df = read_blast_report(report_path)
            qseq = get_top_qseq(df=df, locus=locus)
            if type(qseq) == str:
                locus_length = len(qseq)
            else:
                qseq = "N" * locus_length
            out.write(f"{qseq}\n")


def call_muscle(infile: Path):
    outfile = infile.with_suffix(".align.fasta")
    cmd = f"muscle -in {infile} -out {outfile} -maxiters 1"
    logging.debug(cmd)
    run_subprocess(cmd, get_stdout=True)


def call_fasconcat(target_dir: Path, fasconcat_exec: Path):
    cmd = f"perl {fasconcat_exec} -s -p"
    p = Popen(cmd, shell=True, cwd=str(target_dir))
    p.wait()


if __name__ == "__main__":
    main()
