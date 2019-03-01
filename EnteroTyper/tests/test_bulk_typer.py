import tempfile
import shutil
from EnteroTyper.bin.bulk_typer import get_sample_name_dict, bulk_sample_typing
from pathlib import Path


def test_get_sample_name_dict():
    # Create mock directory structure
    indir = Path(tempfile.mkdtemp())

    fasta_file_1_fasta = tempfile.NamedTemporaryFile(mode='w+', dir=indir, encoding='utf-8', suffix=".fasta")
    fasta_file_1_fasta.write(">Header\nATCG")
    fasta_file_1_fasta.flush()

    fasta_file_2_fna = tempfile.NamedTemporaryFile(mode='w+', dir=indir, encoding='utf-8', suffix=".fna")
    fasta_file_2_fna.write(">Header\nATCG")
    fasta_file_2_fna.flush()

    fasta_file_3_fa = tempfile.NamedTemporaryFile(mode='w+', dir=indir, encoding='utf-8', suffix=".fa")
    fasta_file_3_fa.write(">Header\nATCG")
    fasta_file_3_fa.flush()

    sample_name_dict = get_sample_name_dict(indir=indir)
    assert len(sample_name_dict) == 3

    shutil.rmtree(indir)


def test_bulk_sample_typing():
    test_data_dir = Path(__file__).parent / 'data'
    assert test_data_dir.exists()
    test_database = Path(__file__).parent / 'data' / 'database'
    assert test_database.exists()

    # Set up temporary directory
    tmpoutdir = Path(f'/tmp/typer_testing_bulk')
    if tmpoutdir.exists():
        shutil.rmtree(tmpoutdir)

    detailed_report_list = bulk_sample_typing(indir=test_data_dir, outdir=tmpoutdir, database=test_database,
                                              keep_blast=False)
    # There should be 3 reports generated - one for each input sample
    assert len(detailed_report_list) == 3
