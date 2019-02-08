import click
import logging
import pandas as pd
from pathlib import Path
from enterobase_typer import create_outdir


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    value = Path(value)
    return value


@click.command(help="Takes a list of target *.cgMLST_Allele_Report.tsv files, followed by several options. "
                    "Number of mismatches will be indicated in output files - a zero means no mismatches between types."
               )
@click.option('-o', '--out_dir',
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
def main(targets, out_dir, verbose):
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
    logging.info("Comparing provided sequences")
    call_sequence_comparison(targets=targets, out_dir=out_dir)
    logging.info("Script complete")


def call_sequence_comparison(targets: list, out_dir: Path):
    create_outdir(out_dir=out_dir)

    # Populate target_dict (key=sample name, value=path to report)
    targets = [Path(target) for target in targets]
    target_dict = {}
    for target in targets:
        target_name = target.with_suffix("").name.replace("_cgMLST_Allele_Report", "")
        target_dict[target_name] = target

    # Make pairwise comparisons
    pairwise_dict = {}
    mismatch_dict = {}

    with click.progressbar(target_dict.items(), length=len(target_dict), label='Pairwise comparisons') as bar:
        for target_1, target_1_report in bar:
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

    with click.progressbar(pairwise_dict.items(), length=len(pairwise_dict), label='Writing .xlsx') as bar:
        for key, val in bar:
            df = pd.DataFrame.from_dict(val)
            out_name = out_dir / (key + "_pairwise_comparisons.xlsx")
            if out_name.exists():
                continue
            writer = pd.ExcelWriter(str(out_name), engine='xlsxwriter')
            df.to_excel(writer, sheet_name=key)
            worksheet = writer.sheets[key]
            worksheet.conditional_format('B2:AZ4000', {'type': '2_color_scale'})
            writer.save()

    master_df = []
    out_name = out_dir / "all_samples_comparison.xlsx"
    for key, val in mismatch_dict.items():
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

    logging.debug(f"NUM. MISMATCHES: {mismatches}/{len(difference_dict)}")
    mismatch_pct = (mismatches / len(difference_dict)) * 100
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


if __name__ == "__main__":
    main()
