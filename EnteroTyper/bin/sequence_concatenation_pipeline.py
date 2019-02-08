import os
import logging
import pandas as pd
import multiprocessing
from copy import deepcopy
from tqdm import tqdm
from pathlib import Path
from EnteroTyper.bin.sequence_type_comparison import call_sequence_comparison
from EnteroTyper.bin.enterobase_typer import get_database_files, create_outdir, generate_cgmlst_report
from EnteroTyper.bin.accessories import run_subprocess


def sequence_concatenation_pipeline(targets: list, database: Path, outdir: Path):
    targets = [Path(target) for target in targets]
    database_files = get_database_files(database=database)
    create_outdir(out_dir=outdir)

    report_dict = {}
    for i, target in enumerate(targets):
        target_name = target.with_suffix("").name + f'_{i}'
        report_dict[target_name] = target

    logging.info("Generating FASTA files from reports")
    df_params = [(database_file, outdir, report_dict) for database_file in database_files]
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as p:
        p.starmap(write_fasta, df_params)

    logging.info("Aligning FASTA files with MUSCLE")
    fasta_files = list(outdir.glob("*.fasta"))
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as p:
        p.map(call_muscle, fasta_files)

    # Grab newly aligned fasta files
    aligned_fasta_files = list(outdir.glob("*.align.fasta"))

    # Create multifasta with fasconcat
    logging.info("Concatenating sequences")
    n_processes = int(multiprocessing.cpu_count() / 8)
    sample_ids = list(report_dict.keys())
    outfile = concatenate_sequence_directory(sample_ids=sample_ids, sequence_directory=outdir, n_processes=n_processes,
                                             outdir=outdir)
    logging.info(f"Concatenated sequences available at {outfile}")

    # Remove remaining fasta files
    [os.remove(str(fasta)) for fasta in aligned_fasta_files]

    # Sequence type comparison
    logging.info("Conducting sequence type comparisons")
    sequence_type_report_list = []
    for target in tqdm(targets):
        df = pd.read_csv(target, sep='\t')
        cgmlst_allele_report = generate_cgmlst_report(df=df, out_dir=outdir,
                                                      sample_name=target.name.rsplit("_", 1)[0])
        sequence_type_report_list.append(cgmlst_allele_report)

    call_sequence_comparison(targets=sequence_type_report_list, out_dir=outdir / 'sequence_type_comparisons')

    # Cleanup reports
    for report in sequence_type_report_list:
        try:
            os.remove(str(report))
        except FileNotFoundError:
            continue


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
    """
    Produces an aligned version of an input FASTA file, then removes the original.
    """
    outfile = infile.with_suffix(".align.fasta")
    cmd = f"muscle -in {infile} -out {outfile} -maxiters 1"
    run_subprocess(cmd, get_stdout=True)
    infile.unlink()


def concatenate_sequence_directory(sample_ids: [str], sequence_directory: Path, n_processes: int, outdir: Path) -> Path:
    """
    Given a sequence directory containing aligned multi-FASTA files, will attempt to concatenate all of the sequences
    into a single file. Minimal example showing input and output (concatenated_sequences.fasta):

    files in directory:
        gene_1.fasta
        gene_2.fasta

    contents of gene1.fasta:
        >sample_1
        ATCG
        >sample2
        ATTT

    contents of gene2.fasta:
        >sample_1
        GGGGAGGGGGTCA

    concatenated_sequences.fasta:
        >sample_1
        ATCGGGGGAGGGGGTCA
        >sample_2
        ATTTNNNNNNNNNNNNN

    - Provided sample_ids list must contain exact matches to headers in FASTA files
    - All sequences in a single FASTA file must be the exact same length; this is done through an aligner like MUSCLE
    - Expectation is that every FASTA file will have the same headers, or a subset of headers present in sample_ids
    """
    concat_seqs_dict = generate_concat_seqs_dict(sample_ids=sample_ids, indir=sequence_directory,
                                                 n_processes=n_processes)
    outfile = write_concat_seqs_dict(concat_seqs_dict=concat_seqs_dict, outdir=outdir)
    return outfile


def write_concat_seqs_dict(concat_seqs_dict: dict, outdir: Path) -> Path:
    outfile = outdir / 'concatenated_sequences.fasta'
    if outfile.exists():
        outfile.unlink()
    outfile_ = open(str(outfile), 'a+')
    for sample_id, sequence in concat_seqs_dict.items():
        outfile_.write(f">{sample_id}\n")
        outfile_.write(f"{sequence}\n")
    outfile_.close()
    return outfile


def generate_concat_seqs_dict(sample_ids: set, indir: Path, n_processes=4) -> dict:
    # Potentially takes a lot of RAM. Stores the concatenated sequences for each sample.
    sequence_storage = {sample_id: "" for sample_id in sample_ids}
    fasta_files = sorted(list(indir.glob("*.fasta")))

    # Set # of concurrent processes to run
    pool = multiprocessing.Pool(processes=n_processes)
    cluster_dicts = [
        pool.apply_async(populate_template_dict, args=(sequence_storage, f, sample_ids)) for f in fasta_files
    ]
    cluster_dicts = [result.get() for result in cluster_dicts]

    # Merge all of the dictionaries into one
    for d in cluster_dicts:
        for cluster, sequence in d.items():
            sequence_storage[cluster] += sequence
    return sequence_storage


def populate_template_dict(template_dict: dict, cluster_file: Path, sample_ids: list):
    with open(str(cluster_file), 'r') as f:
        cluster_dict = deepcopy(template_dict)
        lines = f.readlines()
        cluster_samples = []
        seq_length = 0
        seq_lengths = []
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # Add the sequence length to the ongoing tracking list
                if seq_length != 0:
                    seq_lengths.append(seq_length)
                seq_length = 0
                # Remove '>'
                header = line.replace(">", "")
                # Add this header to cluster_samples to keep track of which samples we have sequence for
                cluster_samples.append(header)
            else:
                cluster_dict[header] += line
                seq_length += len(line)
        try:
            assert len(set(seq_lengths)) == 1
        except AssertionError:
            print(f"ERROR: Varying sequence lengths detected in {f.name}")
            print(seq_lengths)
            quit()
        cluster_samples = set(cluster_samples)
        for s in sample_ids:
            if s not in cluster_samples:
                cluster_dict[s] += ('N' * seq_lengths[0])
        return cluster_dict
