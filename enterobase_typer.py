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

    type_sample(input_assembly, database, out_dir, create_db)


def type_sample(input_assembly: Path, database: Path, out_dir: Path, create_db: bool):
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

    # Get reverse complement of minus strands
    df['qseq_strand_aware'] = df.apply(get_reverse_complement_row, axis=1)

    # Drop extraneous columns
    df = df.drop(['qseq', 'sstrand', 'score'], axis=1)

    # Sort
    df = df.sort_values(by=['locus'])

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

    # TODO: Test this blastn move bit of code
    # Move BLASTn files
    blastn_folder = out_dir / 'blastn_output'
    os.makedirs(str(blastn_folder), exist_ok=True)
    blastn_files = list(blastn_folder.glob("*.BLASTn"))
    for f in blastn_files:
        shutil.move(str(f), str(blastn_folder / f.name))

    logging.info(f"=============TYPING COMPLETE=============")
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
    logging.debug(f"Checking {locus_name}...")
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


def get_reverse_complement_row(row):
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


def run_subprocess(cmd: str, get_stdout=False):
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


if __name__ == "__main__":
    main()
