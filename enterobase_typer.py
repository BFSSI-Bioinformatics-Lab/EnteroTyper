#!/usr/bin/env python3
# TODO: Fix the inconsistent naming between the GitHub project and PyCharm project. Rename to EnterobaseTyper.

__version__ = "0.0.1"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"

import os
import gzip
import click
import shutil
import logging
import multiprocessing
import pandas as pd

from numba import jit
from subprocess import Popen, PIPE
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

    cli_call = True
    type_sample(input_assembly=input_assembly, database=database,
                out_dir=out_dir, create_db=create_db, cli_call=cli_call)


def type_sample(input_assembly: Path, database: Path, out_dir: Path, create_db: bool, cli_call: bool = False):
    """
    Function wrapping all of the functionality of the enterobase_typer script. Can be called directly.
    """

    # Output directory creation/validation
    create_outdir(out_dir=out_dir)

    logging.debug(f"input_assembly: {input_assembly}")
    logging.debug(f"database: {database}")
    logging.debug(f"outdir: {out_dir}")
    logging.debug(f"create_db: {create_db}")

    # Call makeblastdb on each database file if create_db=True
    if create_db:
        database_files = makeblastdb_database(database=database, cli_call=cli_call)
    else:
        database_files = get_database_files(database=database)

    # Query input_assembly against each loci in the cgMLST database
    df = multiprocess_blastn_call(database_files, input_assembly, out_dir)

    # Grab sample name from input_assembly
    sample_name = input_assembly.with_suffix("").name.replace(".pilon", "")

    # Prepare detailed report
    detailed_report = generate_detailed_report(df=df, out_dir=out_dir, sample_name=sample_name)

    # Prepare cgMLST report
    cgmlst_allele_report = generate_cgmlst_report(df=df, out_dir=out_dir, sample_name=sample_name)

    # Move BLASTn files
    move_blastn_files(out_dir=out_dir)

    logging.info(f"=============TYPING COMPLETE=============")
    logging.info(f"cgMLST Allele Report: {cgmlst_allele_report}")
    logging.info(f"Detailed Report: {detailed_report}")


def move_blastn_files(out_dir: Path):
    """
    Moves all of the BLASTn output to a new folder 'blastn_output'
    """
    blastn_folder = out_dir / 'blastn_output'
    os.makedirs(str(blastn_folder), exist_ok=True)
    blastn_files = list(out_dir.glob("*.BLASTn"))
    for f in blastn_files:
        shutil.move(str(f), str(blastn_folder / f.name))


def generate_cgmlst_report(df: pd.DataFrame, out_dir: Path, sample_name: str) -> Path:
    """
    Generates the cgMLST report along with a transposed version within the out_dir folder
    """
    cgmlst_allele_report = out_dir / f"{sample_name}_cgMLST_Allele_Report.tsv"
    cgmlst_df = get_sequence_type(df)
    cgmlst_df_transposed = cgmlst_df.transpose()
    cgmlst_df_transposed.to_csv(cgmlst_allele_report, sep="\t", header=False)
    return cgmlst_allele_report


def generate_detailed_report(df: pd.DataFrame, out_dir: Path, sample_name: str) -> Path:
    """
    Generates the detailed report within the out_dir folder
    """
    output_detailed_report = out_dir / f"{sample_name}_BLASTn_Detailed_Report.tsv"
    df.to_csv(output_detailed_report, sep="\t", index=None)
    return output_detailed_report


def get_database_files(database: Path, loci_suffix: str = "*.gz"):
    """
    Grabs database files within a provided directory and returns a list of everything present.
    Matches against a given suffix (defaults to *.gz, the expected database file extension)
    """
    database_files = list(database.glob(loci_suffix))
    if len(database_files) > 0:
        return database_files
    else:
        raise EmptyDatabase(f"Could not find database files at {database}. "
                            f"Try re-running the script with the --create_db flag.")


def create_outdir(out_dir: Path) -> Path:
    """
    Creates output directory, quits/raises error if it already exists
    """
    os.makedirs(str(out_dir), exist_ok=False)
    logging.debug(f"Created directory {out_dir}")
    return out_dir


def multiprocess_blastn_call(database_files: list, input_assembly: Path, outdir: Path) -> [pd.DataFrame]:
    """
    Calls blastn on every database file against the input_assembly, reads top hit of results into a dataframe per loci,
    then returns a single dataframe of the combined results
    """
    df_params = [(database_file, input_assembly, outdir) for database_file in database_files]

    with Pool(multiprocessing.cpu_count() - 1) as p:
        df_list = p.starmap(closest_allele_df, df_params)
    df = combine_dataframes(df_list)

    # Sort values
    df = df.sort_values(by=['sseqid'])

    # Get reverse complement of minus strands
    df['qseq_strand_aware'] = df.apply(get_reverse_complement_row, axis=1)

    # Drop extraneous columns
    df = df.drop(['qseq', 'sstrand', 'score'], axis=1)

    # Sort again
    df = df.sort_values(by=['locus'])

    return df


def closest_allele_df(database_file: Path, input_assembly: Path, outdir: Path) -> pd.DataFrame:
    """
    Queries the input_assembly against a given database file, returns the top hit as a dataframe
    """
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


def process_blastn_df(df: pd.DataFrame, locus_name: str) -> pd.DataFrame:
    """
`   Sorts and filters the DataFrame according to pident and a newly calculated lratio, evaluates top hits, then
    returns a DataFrame with only the top hit along with a 'hit_type' classification.
    :param df: DataFrame containing parsed *.BLASTn results
    :param locus_name: Name of the locus represented in the *.BLASTn file
    :return: df containing only the top hit from the parsed *.BLASTn
    """
    logging.debug(f"Processing {locus_name}...")
    # Filter by length to slen ratio
    df['lratio'] = df['length'] / df['slen']
    df['locus'] = locus_name

    # Filter junk from DataFrame
    df = df.query("lratio >= 0.90 & pident >= 90.0 & lratio <=1.10")

    # Sort values so the best hits are at the top
    df = df.sort_values(["pident", "lratio"], ascending=False)
    df = df.reset_index(drop=True)

    hit_type = None
    if len(df) == 1:
        # CASE 1: PERFECT SINGLE HIT
        if df.pident[0] == 100.0 and df.lratio[0] == 1.0:
            hit_type = "PERFECT_HIT"
    elif len(df) > 1 and df.pident[0] == 100.0 and df.lratio[0] == 1.0:
        # CASE 2: PERFECT DUPLICATE - top two hits both have 100% pident and lratio
        if df.pident[1] == 100.0 and df.lratio[1] == 1.0:
            hit_type = "PERFECT_DUPLICATE"
        # CASE 3: PERFECT HIT w/very close second hit
        elif (0.99 <= df.lratio[1] <= 1.01) or df.pident[1] >= 0.99:
            hit_type = "PERFECT_HIT_WITH_POTENTIAL_PARALOG"
        # CASE 4: PERFECT HIT w/slightly less close second hit
        else:
            hit_type = "PERFECT_HIT"
    elif len(df) == 0:  # no close hits
        hit_type = "NO_MATCH"
        data = []
        data.insert(0, {"qseqid": "NA", "sseqid": "NA", "slen": "NA",
                        "length": "NA", "qstart": "NA", "qend": "NA",
                        "pident": "NA", "score": "NA", "locus": locus_name,
                        "sstrand": "NA", "qseq": "NA", "lratio": "NA"})
        df = pd.concat([pd.DataFrame(data), df], ignore_index=True, sort=False)
    else:
        hit_type = "CLOSEST_HIT"

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
        run_subprocess(cmd, get_stdout=True)
    elif db_file.suffix == ".fasta":
        cmd = f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd, get_stdout=True)
    elif db_file.suffix == "":
        os.rename(str(db_file), str(db_file.with_suffix(".fasta")))
        cmd = f"makeblastdb -in {db_file} -parse_seqids -dbtype nucl -out {db_name} -title {db_name}"
        run_subprocess(cmd, get_stdout=True)
    else:
        logging.debug("Invalid file format provided to call_makeblastdb()")


def makeblastdb_database(database: Path, cli_call: bool, loci_suffix: str = "*.gz") -> list:
    """
    Calls makeblastdb on every *.gz file in a given directory. Intended to function with the files retrieved via the
    EnterobasePull script (https://github.com/bfssi-forest-dussault/EnterobasePull).
    :param database: Path to database retrieved with EnterobasePull
    :param loci_suffix: Suffix of files to run makeblastdb on
    :param cli_call: Flag to show progress bar if this was called via command-line
    """
    logging.debug(f"Calling makeblastdb on database at {database}")
    db_files = list(database.glob(loci_suffix))
    if cli_call:
        with click.progressbar(db_files, length=len(db_files), label='makeblastdb') as bar:
            for f in bar:
                call_makeblastdb(f)
    else:
        for f in db_files:
            call_makeblastdb(f)
    return db_files


def get_allele_length(database_file: Path) -> int:
    with gzip.open(database_file, "rt") as handle:
        lines = [line.strip() for line in handle]
        for line in lines:
            if line[0] != ">":
                return len(line)


def call_blastn(database_file: Path, query_fasta: Path, out_dir: Path) -> (Path, str):
    """
    Calls blastn against a query *.fasta file using the provided Enterobase DB target, dumps output into out_dir
    :param database_file: Formatted (makeblastdb) BLAST database file
    :param query_fasta: Sequence to query against database file
    :param out_dir: Path to directory to store output
    :return: BLASTn output (Path) and the base name of the target being searched (str)
    """
    # Dynamically set word_size
    allele_length = get_allele_length(database_file=database_file.with_suffix(".gz"))
    word_size = round(allele_length / 4)

    locus_name = database_file.with_suffix('').name
    reference_name = query_fasta.with_suffix('').name
    out_file = out_dir / Path(reference_name + "." + locus_name + ".BLASTn")
    cmd = f"blastn -query {query_fasta} -db {database_file} -out {out_file} " \
          f"-outfmt '6 qseqid sseqid slen length qstart qend pident score sstrand qseq' -word_size {word_size}"
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


def get_reverse_complement_row(row) -> str:
    """
    Takes DataFrame row and returns the reverse complement of the sequence if the strand is 'minus'
    :param row:
    :return:
    """
    complement_dict = {'A': 'T',
                       'C': 'G',
                       'G': 'C',
                       'T': 'A'}
    sequence = row['qseq'].upper()
    if row['sstrand'] == 'minus':
        reverse_complement = "".join(complement_dict.get(base, base) for base in reversed(sequence))
        return reverse_complement
    else:
        return sequence


def run_subprocess(cmd: str, get_stdout: bool = False) -> str:
    if get_stdout:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        out = out.decode().strip()
        err = err.decode().strip()
        if out != "":
            return out
        elif err != "":
            return err
        else:
            return ""
    else:
        p = Popen(cmd, shell=True)
        p.wait()


class EmptyDatabase(Exception):
    pass


if __name__ == "__main__":
    main()
