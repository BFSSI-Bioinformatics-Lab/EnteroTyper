#!/usr/bin/env python3
# TODO: Fix the inconsistent naming between the GitHub project and PyCharm project. Rename to EnterbaseTyper.

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
    blastn_file, locus_name = call_blastn(database_file=database_file, query_fasta=input_assembly, outdir=outdir)
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
    # Filter by length to slen ratio
    df['lratio'] = df['length'] / df['slen']
    df['locus'] = locus_name
    df = df.query("lratio >= 0.92 & pident >= 92")  # Filter out anything below 92% lratio and < 92 pident
    df = df.sort_values(by=['pident'], ascending=False)
    df = df.reset_index(drop=True)

    hit_type = "NA"
    if len(df) == 1:
        if df['pident'][0] == 100:  # 100% hit
            hit_type = "PERFECT_MATCH"
            logging.debug(f"PERFECT MATCH\t\tALLELE:{df['sseqid'][0]}\t\tPIDENT:{df['pident'][0]}")
        else:  # high lratio and pident hit, but not in the cgMLST database
            hit_type = "NEW_ALLELE"
            logging.debug(f"NEW ALLELE\t\t\tCLOSEST ALLELE:{df['sseqid'][0]}\t\tPIDENT:{df['pident'][0]}")
    elif len(df) > 1:
        if df.pident[0] == 100 and df.pident[1] == 100:  # weird duplicate 100% hit
            hit_type = "PERFECT_DUPLICATE"
            logging.debug(f"PERFECT DUPLICATE\t\t\tALLELE:{df['sseqid'][0]}\t\tALLELE:{df['sseqid'][1]}")
            df = df.head(1)
        elif df.pident[0] == 100 and df.pident[1] != 100:  # perfect hit + close hit
            hit_type = "LIKELY_PARALOG"
            logging.debug(f"LIKELY PARALOG\t\tPERFECT ALLELE:{df['sseqid'][0]}\t\tCLOSE ALLELE:{df['sseqid'][1]}")
            df = df.head(1)
    elif len(df) == 0:  # no close hit
        hit_type = "NO_MATCH"
        data = []
        data.insert(0, {"qseqid": "NA", "sseqid": "NA", "slen": "NA",
                        "length": "NA", "qstart": "NA", "qend": "NA",
                        "pident": "NA", "score": "NA", "locus": locus_name})
        df = pd.concat([pd.DataFrame(data), df], ignore_index=True, sort=False)
        df = df.head(1)
        logging.debug(f"NO MATCHES WITH LRATIO > 0.92 FOR {locus_name}")

    df["hit_type"] = hit_type
    return df


def call_makeblastdb(db_file: Path):
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
    db_files = list(database.glob(loci_suffix))
    for f in db_files:
        call_makeblastdb(f)


def call_blastn(database_file: Path, query_fasta: Path, outdir: Path) -> tuple:
    locus_name = database_file.with_suffix('').name
    reference_name = query_fasta.with_suffix('').name
    outfile = outdir / Path(reference_name + "." + locus_name + ".BLASTn")
    cmd = f"blastn -query {query_fasta} -db {database_file} -out {outfile} -max_target_seqs 5 " \
          f"-outfmt '6 qseqid sseqid slen length qstart qend pident score sstrand qseq' -word_size 10 -perc_identity 90"
    run_subprocess(cmd)
    return outfile, locus_name


def combine_dataframes(dfs: list) -> pd.DataFrame:
    df = pd.concat(dfs, sort=False)
    return df


def parse_blastn(blastn_file: Path):
    headers = ["qseqid", "sseqid", "slen", "length", "qstart", "qend", "pident", "score", "sstrand", "qseq"]
    df = pd.read_csv(blastn_file, delimiter="\t", names=headers)
    df['lratio'] = df['length'] / df['slen']
    df = df.sort_values(["pident", "lratio"], ascending=False)
    df = df.head(1)  # Only take the top hit
    return df


@jit()
def get_reverse_complement(sequence: str):
    complement_dict = {'A': 'T',
                       'C': 'G',
                       'G': 'C',
                       'T': 'A'}
    reverse_complement = "".join(complement_dict.get(base, base) for base in reversed(sequence))
    return reverse_complement


def run_subprocess(cmd: str):
    p = Popen(cmd, shell=True)
    p.wait()


if __name__ == "__main__":
    main()
