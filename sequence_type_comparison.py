from pathlib import Path
import pandas as pd


def compare_reports(report_1: Path, report_2: Path) -> dict:
    dict_1 = read_allele_report(report_1)
    dict_2 = read_allele_report(report_2)
    difference_dict = compare_allele_dicts(dict_1, dict_2)
    mismatches = count_mismatches(difference_dict)
    # print("DIFFERENCE DICTIONARY:")
    # for key, val in difference_dict.items():
    #     print(f"{key}: {val}")
    print(f"NUM. MISMATCHES: {mismatches}/{len(difference_dict)}")
    mismatch_pct = (mismatches/len(difference_dict)) * 100
    # print(f"PERCENT MISMATCH: {mismatch_pct}%")
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
        print(f"WARNING: Length of dictionaries are not the same. {len(keys_1)} vs. {len(keys_2)}")

    for key in keys_1:
        if key not in keys_2:
            error_state = True
            print(f"WARNING: {key} not in both dictionaries")

    for key in keys_2:
        if key not in keys_1:
            error_state = True
            print(f"WARNING: {key} not in both dictionaries")

    if error_state:
        print(f"FAIL: Provided dictionaries are not equivalent. Beware.")
    else:
        pass
        # print(f"SUCCESS: Dictionaries are equivalent")


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
