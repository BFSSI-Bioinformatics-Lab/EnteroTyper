import tempfile

from EnteroTyper.bin.typer import *
from pathlib import Path


def test_type_sample():
    """
    Uses dummy data in the test_data_dir directory. Uses the MLST_Achtman database.
    assembly_1.fasta should contain only recA_1 and fumC_7
    assembly_2.fasta should contain only recA_1 and fumC_16
    assembly_3.fasta should not contain any markers
    :return:
    """
    # Establish test data directories, make sure they are structured correctly
    test_data_dir = Path(__file__).parent / 'data'
    assert test_data_dir.exists()
    test_database = Path(__file__).parent / 'data' / 'database'
    assert test_database.exists()

    test_assembly_1 = test_data_dir / 'assembly_1.fasta'  # recA_1 and fumC_7
    test_assembly_2 = test_data_dir / 'assembly_2.fasta'  # recA_1 and fumC_16
    test_assembly_3 = test_data_dir / 'assembly_3.fasta'  # No markers
    test_assemblies = [test_assembly_1, test_assembly_2, test_assembly_3]

    # Test if reports are generated
    for assembly in test_assemblies:
        assert assembly.exists()
        # Set up temporary directory
        tmpoutdir = Path(f'/tmp/typer_testing_{assembly.name}')
        if tmpoutdir.exists():
            shutil.rmtree(tmpoutdir)
        # Run typer
        detailed_report = type_sample(input_assembly=assembly, database=test_database,
                                      outdir=tmpoutdir, create_db=True, keep_blast=False)
        # Ensure report was generated
        assert detailed_report.exists()

        # Verify contents of report
        df = pd.read_csv(detailed_report, sep="\t")
        expected_columns = ['qseqid', 'sseqid', 'slen', 'length', 'qstart', 'qend',
                            'pident', 'lratio', 'locus', 'hit_type', 'qseq_strand_aware']
        if 'assembly_3' not in assembly.name:
            assert list(df.columns) == expected_columns

        # Verify individual assembly test case results
        if 'assembly_1' in assembly.name:
            assert 'fumC_7' in df['sseqid'].tolist()
            assert 'recA_1' in df['sseqid'].tolist()
        elif 'assembly_2' in assembly.name:
            assert 'fumC_16' in df['sseqid'].tolist()
            assert 'recA_1' in df['sseqid'].tolist()
        elif 'assembly_3' in assembly.name:
            hit_types = list(df['hit_type'].unique())
            assert len(hit_types) == 1
            assert hit_types[0] == 'NO_MATCH'

        # Remove temporary directory
        shutil.rmtree(tmpoutdir)


def test_move_blastn_files():
    # Create mock directory structure
    outdir = Path(tempfile.mkdtemp())
    indir = Path(tempfile.mkdtemp())
    blastn_file = tempfile.NamedTemporaryFile(mode='w+', dir=indir, encoding='utf-8', suffix=".BLASTn")
    blastn_file.write("tmp")
    blastn_file.flush()

    # Ensure indir has BLASTn file inside of it
    assert len((list(indir.glob("*.BLASTn")))) != 0

    # Call method on indir and confirm that BLASTn file was moved to outdir
    move_blastn_files(indir=indir, outdir=outdir)
    assert len((list(outdir.glob("*.BLASTn")))) != 0

    # Cleanup
    shutil.rmtree(outdir)
    shutil.rmtree(indir)


def test_delete_blastn_files():
    # Create mock directory structure
    indir = Path(tempfile.mkdtemp())
    blastn_file = tempfile.NamedTemporaryFile(mode='w+', dir=indir, encoding='utf-8', suffix=".BLASTn")
    blastn_file.write("tmp")
    blastn_file.flush()

    # Ensure indir has BLASTn file inside of it
    assert len((list(indir.glob("*.BLASTn")))) != 0
    delete_blastn_files(indir=indir)

    # Ensure file was deleted
    assert len((list(indir.glob("*.BLASTn")))) == 0


def test_generate_cgmlst_report():
    pass


def test_generate_detailed_report():
    pass


def test_get_database_files():
    pass


def test_create_outdir():
    pass


def test_multiprocess_blastn_call():
    pass


def test_closest_allele_df():
    pass


def test_get_sequence_type():
    pass


def test_process_blastn_df():
    pass


def test_call_makeblastdb():
    pass


def test_makeblastdb_database():
    pass


def test_get_allele_length():
    pass


def test_call_blastn():
    pass


def test_combine_dataframes():
    pass


def test_parse_blastn():
    pass


def test_get_reverse_complement_row():
    pass
