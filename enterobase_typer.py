#!/usr/bin/env python3
# TODO: Fix the inconsistent naming between the GitHub project and PyCharm project. Rename to EnterobaseTyper.

__version__ = "0.0.1"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"

import os
import click
import logging
import multiprocessing
import pandas as pd

from numba import jit
from subprocess import Popen
from multiprocessing import Pool

from pathlib import Path

script = os.path.basename(__file__)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    logging.info(f"Version: {__version__}")
    logging.info(f"Author: {__author__}")
    logging.info(f"Email: {__email__}")
    quit()


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    value = Path(value)
    return value


@click.command()
@click.option('-i', '--input_assembly',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Path to input assembly in FASTA format',
              callback=convert_to_path)
@click.option('-db', '--database',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Path to your MLST database',
              callback=convert_to_path)
@click.option('-o', '--out_dir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Root directory to store all output files',
              callback=convert_to_path)
@click.option('--create_db',
              help='Set this flag to create the blastDB files using makeblastdb in the specified database directory.'
                   'Will re-create the database files if they are already present.',
              is_flag=True,
              required=False,
              default=False)
@click.option('-v', '--verbose',
              is_flag=True,
              default=False,  # Set this to false eventually
              help='Set this flag to enable more verbose logging.')
@click.option('--version',
              help='Specify this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def main(input_assembly: Path, database: Path, out_dir: Path, create_db: bool, verbose: bool):
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

    # Output directory validation
    try:
        os.makedirs(str(out_dir), exist_ok=False)
        logging.debug(f"Created directory {out_dir}")
    except FileExistsError:
        logging.error("ERROR: Output directory already exists.")
        quit()

    logging.debug(f"input_assembly: {input_assembly}")
    logging.debug(f"database: {database}")
    logging.debug(f"outdir: {out_dir}")
    logging.debug(f"create_db: {create_db}")

    if create_db:
        logging.debug(f"Calling makeblastdb on database at {database}")
        makeblastdb_database(database=database)
    database_files = list(database.glob("*.gz"))

    df_list = prepare_df_list_multiprocess(database_files, input_assembly, out_dir)
    df = combine_dataframes(df_list)
    df = df.sort_values(by=['sseqid'])

    # Prepare detailed report
    output_detailed_report = out_dir / "BLASTn_Detailed_Report.tsv"
    df.to_csv(output_detailed_report, sep="\t", index=None)

    # Prepare cgMLST report
    cgmlst_allele_report = out_dir / "cgMLST_Allele_Report.tsv"
    cgmlst_allele_report_transposed = out_dir / "cgMLST_Allele_Report_transposed.tsv"
    cgmlst_df = get_sequence_type(df)
    cgmlst_df_transposed = cgmlst_df.transpose()
    cgmlst_df.to_csv(cgmlst_allele_report, sep="\t", index=None)
    cgmlst_df_transposed.to_csv(cgmlst_allele_report_transposed, sep="\t", header=False)

    logging.info(f"Done!")
    logging.info(f"cgMLST Allele Report: {cgmlst_allele_report}")
    logging.info(f"Detailed Report: {output_detailed_report}")


def prepare_df_list_multiprocess(database_files: list, input_assembly: Path, outdir: Path):
    df_params = [(database_file, input_assembly, outdir) for database_file in database_files]
    with Pool(multiprocessing.cpu_count() - 1) as p:
        df_list = p.starmap(closest_allele_df, df_params)
    return df_list


def closest_allele_df(database_file, input_assembly, outdir):
    database_file = database_file.with_suffix(".blastDB")
    blastn_file, locus_name = call_blastn(database_file=database_file, query_fasta=input_assembly, out_dir=outdir)
    df = parse_blastn(blastn_file)
    df = process_blastn_df(df, locus_name)
    return df


@jit()
def get_sequence_type(df: pd.DataFrame) -> pd.DataFrame:
    loci = df['locus'].unique()
    sequence_type_dict = {}
    for locus in loci:
        sequence_type_df = df[df['locus'] == locus]
        sequence_type = sequence_type_df['sseqid'].iloc[0]
        if str(sequence_type) != 'nan':
            sequence_type = sequence_type.rsplit('_', 1)[-1]
        else:
            sequence_type = "NA"
        sequence_type_dict[locus] = sequence_type
    df = pd.DataFrame.from_dict(sequence_type_dict, orient='index').transpose()
    return df


def process_blastn_df(df: pd.DataFrame, locus_name: str):
    """
`   Sorts and filters the DataFrame according to pident and a newly calculated lratio, evaluates top hits, then
    returns a DataFrame with only the top hit along with a 'hit_type' classification.
    :param df: DataFrame containing parsed *.BLASTn results
    :param locus_name: Name of the locus represented in the *.BLASTn file
    :return: df containing only the top hit from the parsed *.BLASTn
    """
    # Filter by length to slen ratio
    df['lratio'] = df['length'] / df['slen']
    df['locus'] = locus_name
    df = df.query("lratio >= 0.99 & pident >= 99.9 & lratio <=1.00")  # Filter DataFrame

    # Sort values so the best hits are at the top
    df = df.sort_values(["pident", "lratio"], ascending=False)
    df = df.reset_index(drop=True)

    hit_type = "NA"
    if len(df) == 1:
        if df.pident[0] == 100.0 and df.lratio[0] == 1.0:  # 100% hit
            hit_type = "PERFECT_SINGLE_HIT"
            logging.debug(f"{hit_type}\t\tALLELE:{df.sseqid[0]}\t\tPIDENT:{df.pident[0]}")
        else:  # high lratio and pident hit, but not in the cgMLST database
            hit_type = "NEW_ALLELE"
            logging.debug(f"{hit_type}\t\t\tCLOSEST ALLELE:{df.sseqid[0]}\t\tPIDENT:{df.pident[0]}")
    elif len(df) > 1:
        if df.pident[0] == 100.0 and df.pident[1] == 100.0 and df.lratio[0] == 1.0 and df.lratio[1] == 1.0:
            hit_type = "PERFECT_DUPLICATE"
            logging.debug(f"{hit_type}\t\t\tALLELE:{df.sseqid[0]}\t\tALLELE:{df.sseqid[1]}")
        # TODO: Think this is mistakenly being called too often as a result of blastn word_size. Debug.
        elif df.pident[0] == 100.0 and df.pident[1] != 100.0 and df.lratio[0] == 1.0:  # perfect hit + close hit
            hit_type = "PERFECT_HIT_WITH_PARALOG"
            logging.debug(
                f"{hit_type}\t\tPERFECT_HIT:{df.sseqid[0]}\t\t"
                f"CLOSEST_ALLELE:{df.sseqid[1]}({df.pident[1]},{df.lratio[1]})")
    elif len(df) == 0:  # no close hit
        hit_type = "NO_MATCH"
        data = []
        data.insert(0, {"qseqid": "NA", "sseqid": "NA", "slen": "NA",
                        "length": "NA", "qstart": "NA", "qend": "NA",
                        "pident": "NA", "score": "NA", "locus": locus_name,
                        "sstrand": "NA", "qseq": "NA", "lratio": "NA"})
        df = pd.concat([pd.DataFrame(data), df], ignore_index=True, sort=False)
        logging.debug(f"{hit_type} FOR {locus_name}")
    df["hit_type"] = hit_type
    return df.head(1)


def call_makeblastdb(db_file: Path):
    """
    Makes a system call to makeblastdb on a given database file. Can handle *.gz, *.fasta, or no suffix.
    :param db_file: Path to database file
    """
    db_name = db_file.with_suffix(".blastDB")
    if db_file.suffix == ".gz":
        cmd = f"gunzip -c {db_file} | makeblastdb -in - -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd)
    elif db_file.suffix == ".fasta":
        cmd = f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd)
    elif db_file.suffix == "":
        os.rename(str(db_file), str(db_file.with_suffix(".fasta")))
        cmd = f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd)
    else:
        logging.debug("Invalid file format provided to call_makeblastdb()")


def makeblastdb_database(database: Path, loci_suffix: str = "*.gz"):
    """
    Calls makeblastdb on every *.gz file in a given directory. Intended to function with the files retrieved via the
    EnterobasePull script (https://github.com/bfssi-forest-dussault/EnterobasePull).
    :param database: Path to database retrieved with EnterobasePull
    :param loci_suffix: Suffix of files to run makeblastdb on
    """
    db_files = list(database.glob(loci_suffix))
    for f in db_files:
        call_makeblastdb(f)


def call_blastn(database_file: Path, query_fasta: Path, out_dir: Path) -> (Path, str):
    """
    Calls blastn against a query *.fasta file using the provided Enterobase DB target, dumps output into out_dir
    :param database_file: Formatted (makeblastdb) BLAST database file
    :param query_fasta: Sequence to query against database file
    :param out_dir: Path to directory to store output
    :return: BLASTn output (Path) and the base name of the target being searched (str)
    """
    locus_name = database_file.with_suffix('').name
    reference_name = query_fasta.with_suffix('').name
    out_file = out_dir / Path(reference_name + "." + locus_name + ".BLASTn")
    cmd = f"blastn -query {query_fasta} -db {database_file} -out {out_file} " \
          f"-outfmt '6 qseqid sseqid slen length qstart qend pident score sstrand qseq' -word_size 10 " \
          f"-perc_identity 97 -dust yes"
    run_subprocess(cmd)
    return out_file, locus_name


def combine_dataframes(dfs: [pd.DataFrame]) -> pd.DataFrame:
    """
    Receives a list of DataFrames and concatenates them. They must all have the same header.
    :param dfs: List of DataFrames
    :return: Single concatenated DataFrame
    """
    df = pd.concat(dfs, sort=False)
    return df


def parse_blastn(blastn_file: Path) -> pd.DataFrame:
    """
    Parses *.BLASTn file generated by call_blastn(), then returns the df
    :param blastn_file: file path to *.BLASTn file
    :return: DataFrame contents of *.BLASTn file
    """
    headers = ["qseqid", "sseqid", "slen", "length", "qstart", "qend", "pident", "score", "sstrand", "qseq"]
    df = pd.read_csv(blastn_file, delimiter="\t", names=headers)
    return df


@jit()
def get_reverse_complement(sequence: str) -> str:
    """
    Unsophisticated function to return the reverse complement of a DNA string
    :param sequence: target sequence to complement
    :return: complemented sequence
    """
    complement_dict = {'A': 'T',
                       'C': 'G',
                       'G': 'C',
                       'T': 'A'}
    sequence = sequence.upper()
    reverse_complement = "".join(complement_dict.get(base, base) for base in reversed(sequence))
    return reverse_complement


def run_subprocess(cmd: str):
    """
    Makes an external system call and runs it via shell
    :param cmd: string containing system command
    """
    p = Popen(cmd, shell=True)
    p.wait()


if __name__ == "__main__":
    main()
