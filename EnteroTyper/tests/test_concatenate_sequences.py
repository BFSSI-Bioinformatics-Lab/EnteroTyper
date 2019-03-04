from pathlib import Path
import multiprocessing
import shutil
import tempfile
import pandas as pd
import EnteroTyper.bin.concatenate_sequences as cs
import EnteroTyper.bin.bulk_typer as bt
import EnteroTyper.bin.typer as ty

# Setup basic test data structure
test_data_dir = Path(__file__).parent / 'data'
test_database = Path(__file__).parent / 'data' / 'database'
assert test_data_dir.exists()
assert test_database.exists()


def test_sequence_concatenation_pipeline():
    outdir = Path(tempfile.mkdtemp())
    tmpoutdir = Path(f'/tmp/sequence_concat_pipeline_test')
    if tmpoutdir.exists():
        shutil.rmtree(tmpoutdir)
    detailed_report_list = bt.bulk_sample_typing(indir=test_data_dir,
                                                 outdir=outdir,
                                                 database=test_database,
                                                 keep_blast=False,
                                                 create_db=True)
    cs.sequence_concatenation_pipeline(targets=detailed_report_list,
                                       database=test_database,
                                       outdir=tmpoutdir)

    # TODO: Think of some better asserts. Obviously the tmpoutdir exists.
    assert tmpoutdir.exists()

    shutil.rmtree(outdir)


def test_get_top_qseq():
    d = {'locus': ['test1', 'test2'], 'qseq_strand_aware': ['ATCG', 'NA']}
    df = pd.DataFrame(data=d)
    assert cs.get_top_qseq(df=df, locus='test1') == 'ATCG'
    assert cs.get_top_qseq(df=df, locus='test2') == ''


def test_write_fasta():
    # Preliminary setup
    outdir = Path(tempfile.mkdtemp())
    detailed_report_list = bt.bulk_sample_typing(indir=test_data_dir,
                                                 outdir=outdir,
                                                 database=test_database,
                                                 keep_blast=False)
    report_dict = cs.build_report_dict(targets=detailed_report_list)

    # Grab a database file
    database_files = ty.get_database_files(database=test_database)
    database_file = database_files[0]

    # Test write_fasta()
    outfile = cs.write_fasta(database_file=database_file, outdir=outdir, report_dict=report_dict)
    assert outfile.exists()

    # Cleanup
    shutil.rmtree(outdir)


def test_build_report_dict():
    outdir = Path(tempfile.mkdtemp())
    detailed_report_list = bt.bulk_sample_typing(indir=test_data_dir,
                                                 outdir=outdir,
                                                 database=test_database,
                                                 keep_blast=False)
    report_dict = cs.build_report_dict(targets=detailed_report_list)
    assert len(report_dict) == 3
    shutil.rmtree(outdir)


def test_call_muscle():
    outdir = Path(tempfile.mkdtemp())

    # Create copy of assembly_1.fasta and run MUSCLE on it
    test_assembly_1 = test_data_dir / 'assembly_1.fasta'  # recA_1 and fumC_7
    infile = outdir / 'assembly_1_copy.fasta'
    shutil.copy(test_assembly_1, infile)
    outfile = cs.call_muscle(infile)

    # Assert the aligned file exists
    assert outfile.exists()

    # Cleanup
    shutil.rmtree(outdir)


# def test_concatenate_sequence_directory():
#     # Must be equivalent to the FASTA headers of the input assemblies
#     sample_ids = [
#         "Assembly1-recA_1", "Assembly1-fumC_7",  # assembly_1
#         "Assembly2-recA_1", "Assembly2-fumC_16",  # assembly_2
#         "Assembly3-Junk"  # assembly_3
#     ]
#     n_processes = multiprocessing.cpu_count()
#     sequence_directory = Path(tempfile.mkdtemp())
#     outdir = Path(tempfile.mkdtemp())
#
#     # Copy all of the test datafiles to a tmpdir, then align them with muscle
#     test_fastas = list(test_data_dir.glob("*.fasta"))
#     [shutil.copy(f, sequence_directory / f.name) for f in test_fastas]
#     [cs.call_muscle(f) for f in list(sequence_directory.glob("*.fasta"))]
#
#     cs.concatenate_sequence_directory(sample_ids=sample_ids,
#                                       sequence_directory=sequence_directory,
#                                       n_processes=n_processes,
#                                       outdir=outdir)
#     assert outdir.exists()
#     shutil.rmtree(sequence_directory)


def test_write_concat_seqs_dict():
    pass


def test_generate_concat_seqs_dict():
    pass


def test_populate_template_dict():
    pass
