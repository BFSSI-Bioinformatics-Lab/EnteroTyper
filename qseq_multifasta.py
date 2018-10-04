import multiprocessing
from pathlib import Path
from enterobase_typer import run_subprocess
import pandas as pd

"""
cd /mnt/QuizBoy/Elton_McGill_6IsolateRequest_September2018/mlst_typing/FASconcat

for x in *.fasta; do muscle -in $x -out $x.align -maxiters 1; done

perl /home/brock/PycharmProjects/cgMLST_Typing/FASconCAT-G_v1.04.pl -s -p
"""


# TODO: Make this a generic script instead of this hardcoded stuff. Will probably be doing this again.


def main():
    targets = [
        '2013AM0055',
        '2013AM1918',
        '2014AM3028',
        '2014AM-2863',
        'FSIS1502169',
        'FSIS1502916',
        'FSIS1502967',
        'FSIS1502973',
        'FSIS1504606',
        'N55391',
        'CVM19633',  # This sample was not assembled via my pipeline despite the filename ".pilon."
        'G3A',
        'G10A',
        'G11A',
        'G12A',
        'G13A',
        'G15A'
    ]

    report_dict = {}
    for target in targets:
        report_dict[target] = Path(
            f"/mnt/QuizBoy/Elton_McGill_6IsolateRequest_September2018/mlst_typing/{target}/BLASTn_Detailed_Report.tsv")

    database = Path("/mnt/Dean-Venture/EnteroBase/cgMLST_v2")
    database_files = list(database.glob("*.gz"))

    outdir = Path("/mnt/QuizBoy/Elton_McGill_6IsolateRequest_September2018/mlst_typing/FASconcat")

    df_params = [(database_file, outdir, report_dict) for database_file in database_files]
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as p:
        p.starmap(write_fasta, df_params)

    fasta_files = list(outdir.glob("*.fasta"))
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as p:
        p.map(call_muscle, fasta_files)


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
    print(f"{locus}...")
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
    run_subprocess(cmd)




if __name__ == "__main__":
    main()
