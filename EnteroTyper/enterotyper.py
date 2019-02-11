import click
import logging
from pathlib import Path
from EnteroTyper.bin.typer import type_sample, makeblastdb_database
from EnteroTyper.bin.bulk_typer import bulk_sample_typing
from EnteroTyper.bin.concatenate_sequences import sequence_concatenation_pipeline
from EnteroTyper.bin.compare_sequences_types import call_sequence_comparison
from EnteroTyper.bin.accessories import check_all_dependencies

from EnteroTyper.__init__ import __email__, __author__, __version__

DEPENDENCIES = [
    'blastn',
    'blastx',
    'makeblastdb',
    'muscle'
]


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    value = Path(value)
    return value


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    logging.info(f"Version: {__version__}")
    logging.info(f"Author: {__author__}")
    logging.info(f"Email: {__email__}")
    quit()


@click.group()
def enterotyper():
    pass


@enterotyper.command(help="Types an assembly based on an Enterobase scheme",
                     short_help="Types an assembly based on an Enterobase scheme")
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
@click.option('-o', '--outdir',
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
def typer(input_assembly: Path, database: Path, outdir: Path, create_db: bool, verbose: bool):
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

    logging.info("Starting Enterotyper typer script")
    check_all_dependencies(DEPENDENCIES)
    type_sample(input_assembly=input_assembly, database=database,
                outdir=outdir, create_db=create_db)
    logging.info("Script complete")


@enterotyper.command(help="Runs typer on many samples simultaneously",
                     short_help="Runs typer on many samples simultaneously")
@click.option('-i', '--indir',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Path to directory containing FASTA assemblies',
              callback=convert_to_path)
@click.option('-db', '--database',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Path to your MLST database',
              callback=convert_to_path)
@click.option('-o', '--outdir',
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
def bulk(indir: Path, database: Path, outdir: Path, create_db: bool, verbose: bool):
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
    logging.info("Starting Enterotyper bulk script")
    check_all_dependencies(DEPENDENCIES)

    # Set up the database if needed
    if create_db:
        makeblastdb_database(database=database)

    # Type every sample
    detailed_report_list = bulk_sample_typing(indir=indir, outdir=outdir, database=database)

    # Concatenate sequences
    sequence_concatenation_pipeline(targets=detailed_report_list,
                                    outdir=outdir / 'concatenated_sequences',
                                    database=database)
    logging.info("Script complete")


@enterotyper.command(short_help="Concatenate sequences from report files generated by typer",
                     help="Takes a list of target *.BLASTn_Detailed_Report.tsv files as input. "
                          "Extracts sequences from the BLASTn report files as FASTA files, aligns them all with "
                          "MUSCLE, and then concatenates all sequences into a single FASTA. ")
@click.option('-o', '--outdir',
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
def concatenate(targets, outdir, database, verbose):
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

    logging.info("Starting Enterotyper concatenate script")
    check_all_dependencies(DEPENDENCIES)
    sequence_concatenation_pipeline(targets=targets, database=database, outdir=outdir)
    logging.info("Script complete")


@enterotyper.command(short_help="Compares sequence types across samples processed by typer",
                     help="Takes a list of target *.cgMLST_Allele_Report.tsv files as input. "
                          "Number of mismatches will be indicated in output files - a zero means no mismatches "
                          "between types.")
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Root directory to store all output files',
              callback=convert_to_path)
@click.option('-v', '--verbose',
              is_flag=True,
              default=False,
              help='Set this flag to enable more verbose logging.')
@click.argument('targets', nargs=-1, type=click.Path(exists=True))
def compare(targets, outdir, verbose):
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
    logging.info("Starting Enterotyper compare script")
    check_all_dependencies(DEPENDENCIES)
    call_sequence_comparison(targets=targets, outdir=outdir)
    logging.info("Script complete")


if __name__ == "__main__":
    cli = click.CommandCollection(sources=[enterotyper],
                                  help="Enterotyper is a suite of tools to type and evaluate assemblies with "
                                       "Enterobase typing schemes.")
    cli()
