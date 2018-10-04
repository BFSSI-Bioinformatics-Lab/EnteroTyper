from sequence_type_comparison import compare_reports, count_mismatches
from pathlib import Path
import pandas as pd

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
        'CVM19633',
        'G3A',
        'G10A',
        'G11A',
        'G12A',
        'G13A',
        'G15A'
    ]
    type_folder = Path("/mnt/QuizBoy/Elton_McGill_6IsolateRequest_September2018/mlst_typing")
    outdir = Path("/mnt/QuizBoy/Elton_McGill_6IsolateRequest_September2018/mlst_typing/pairwise_comparisons")

    # Make pairwise comparisons
    pairwise_dict = {}
    mismatch_dict = {}
    for target_1 in targets:
        target_1_report = type_folder / target_1 / 'cgMLST_Allele_Report_transposed.tsv'
        inner_dict = {}
        inner_mismatch_dict = {}
        for target_2 in targets:
            target_2_report = type_folder / target_2 / 'cgMLST_Allele_Report_transposed.tsv'

            # Difference dict
            difference_dict = compare_reports(target_1_report, target_2_report)
            inner_dict[target_2] = difference_dict

            # Total mismatches
            mismatches = count_mismatches(difference_dict)
            inner_mismatch_dict[target_2] = mismatches

        mismatch_dict[target_1] = inner_mismatch_dict
        pairwise_dict[target_1] = inner_dict

    for key, val in pairwise_dict.items():
        print(key)
        df = pd.DataFrame.from_dict(val)
        outname = outdir / (key + "_pairwise_comparisons.xlsx")
        if outname.exists():
            continue
        writer = pd.ExcelWriter(str(outname), engine='xlsxwriter')
        df.to_excel(writer, sheet_name=key)
        worksheet = writer.sheets[key]
        worksheet.conditional_format('B2:AZ4000', {'type': '2_color_scale'})
        writer.save()

    master_df = []
    outname = outdir / "all_samples_comparison.xlsx"
    for key, val in mismatch_dict.items():
        df = pd.DataFrame.from_dict(val, orient='index').transpose()
        df['id'] = key
        df = df.set_index('id')
        master_df.append(df)
    df = pd.concat(master_df)

    writer = pd.ExcelWriter(str(outname), engine='xlsxwriter')
    df.to_excel(writer, sheet_name='all_samples_comparison')
    worksheet = writer.sheets['all_samples_comparison']
    worksheet.conditional_format('B2:AZ200', {'type': '3_color_scale'})
    writer.save()


if __name__ == "__main__":
    main()
