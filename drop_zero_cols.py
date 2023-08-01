import pandas as pd
import sys


def zero_col_dropper(input_file):

    df = pd.read_csv(input_file, sep="\t", dtype={"name/position": str, "outcome": int, "*": int})

    zero_sum_cols = []

    for col in df.columns:
        if col == "name/position":
            continue
        if col == "outcome":
            continue
        else:
            if df[col].sum() == 0:
                zero_sum_cols.append(col)

    dropped_df = df.drop(zero_sum_cols, axis=1)

    dropped_df.to_csv(input_file[:-4] + "_dropped_zero_cols.tsv", sep="\t", index=False)
