import logging
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from EnteroTyper.bin.enterobase_typer import create_outdir


def call_sequence_comparison(targets: list, outdir: Path):
    create_outdir(outdir=outdir)

    # Populate target_dict (key=sample name, value=path to report)
    targets = [Path(target) for target in targets]
    target_dict = {}
    for target in targets:
        target_name = target.with_suffix("").name.replace("_cgMLST_Allele_Report", "")
        target_dict[target_name] = target

    # Make pairwise comparisons
    pairwise_dict = {}
    mismatch_dict = {}

    logging.info("Counting mismatches between samples")
    for target_1, target_1_report in tqdm(target_dict.items()):
        inner_dict = {}
        inner_mismatch_dict = {}
        for target_2, target_2_report in target_dict.items():
            # Difference dict
            difference_dict = compare_reports(target_1_report, target_2_report)
            inner_dict[target_2] = difference_dict

            # Total mismatches
            mismatches = count_mismatches(difference_dict)
            inner_mismatch_dict[target_2] = mismatches

        mismatch_dict[target_1] = inner_mismatch_dict
        pairwise_dict[target_1] = inner_dict

    logging.info("Conducting pairwise comparisons")
    for key, val in tqdm(pairwise_dict.items()):
        df = pd.DataFrame.from_dict(val)
        out_name = outdir / (key + "_pairwise_comparisons.xlsx")
        if out_name.exists():
            continue
        writer = pd.ExcelWriter(str(out_name), engine='xlsxwriter')
        df.to_excel(writer, sheet_name=key)
        worksheet = writer.sheets[key]
        worksheet.conditional_format('B2:AZ4000', {'type': '2_color_scale'})
        writer.save()

    master_df = []
    out_name = outdir / "all_samples_comparison.xlsx"
    logging.info("Processing data for final comparison report")
    for key, val in tqdm(mismatch_dict.items()):
        df = pd.DataFrame.from_dict(val, orient='index').transpose()
        df['id'] = key
        df = df.set_index('id')
        master_df.append(df)
    df = pd.concat(master_df)

    writer = pd.ExcelWriter(str(out_name), engine='xlsxwriter')
    df.to_excel(writer, sheet_name='all_samples_comparison')
    worksheet = writer.sheets['all_samples_comparison']
    worksheet.conditional_format('B2:AZ200', {'type': '3_color_scale'})
    writer.save()


def compare_reports(report_1: Path, report_2: Path) -> dict:
    dict_1 = read_allele_report(report_1)
    dict_2 = read_allele_report(report_2)
    difference_dict = compare_allele_dicts(dict_1, dict_2)
    mismatches = count_mismatches(difference_dict)
    mismatch_pct = (mismatches / len(difference_dict)) * 100
    logging.debug(f"NUM. MISMATCHES: {mismatches}/{len(difference_dict)}")
    logging.debug(f"PERCENT MISMATCH: {mismatch_pct}%")
    return difference_dict


def read_allele_report(report: Path) -> dict:
    df = pd.read_csv(report, delimiter="\t", index_col=None, header=None)
    dct = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
    for key, val in dct.items():
        try:
            dct[key] = str(int(val))
        except ValueError:
            dct[key] = "NA"
    return dct


def compare_keys(dict_1: dict, dict_2: dict):
    keys_1 = dict_1.keys()
    keys_2 = dict_2.keys()
    error_state = False

    if len(keys_1) != len(keys_2):
        error_state = True
        logging.warning(f"WARNING: Length of dictionaries are not the same. {len(keys_1)} vs. {len(keys_2)}")

    for key in keys_1:
        if key not in keys_2:
            error_state = True
            logging.warning(f"WARNING: {key} not in both dictionaries")

    for key in keys_2:
        if key not in keys_1:
            error_state = True
            logging.warning(f"WARNING: {key} not in both dictionaries")

    if error_state:
        logging.error(f"FAIL: Provided dictionaries are not equivalent. Beware.")
    else:
        pass


def compare_allele_dicts(dict_1: dict, dict_2: dict):
    # Error check
    compare_keys(dict_1, dict_2)

    # Produce match dictionary
    difference_dict = {}
    for key, val in dict_1.items():
        if key not in dict_2.keys():
            difference_dict[key] = 1  # Missing (counts as no match)
            continue

        if dict_1[key] == dict_2[key]:
            difference_dict[key] = 0  # Match
        else:
            difference_dict[key] = 1  # No match
    return difference_dict


def count_mismatches(difference_dict: dict) -> int:
    mismatch_counter = 0
    for key, val in difference_dict.items():
        mismatch_counter += val
    return mismatch_counter
