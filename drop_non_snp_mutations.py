import pandas as pd

def drop_non_snp_columns(input_file, antibiotics):

    df = pd.read_csv(input_file, sep="\t",
                     dtype={"name/position": str, "outcome": str})
    
    cols_to_drop = []

    for col in df.columns:
        if col in antibiotics:
            continue
        elif col == "name/position":
            continue
        elif col == "outcome":
            continue
        else:
            try:
                temp_col = col[1:-1]
                split_col = [x.strip() for x in temp_col.split(',')]
                split_col[1] = split_col[1][1:]
                split_col[2] = split_col[2][:-1]
                split_col[3] = split_col[3][1:-1]
                type = split_col[3]
                if type == "snp":
                    continue
                else:
                    cols_to_drop.append(col)
            except:
                continue

    df = df.drop(cols_to_drop, axis=1)

    df.to_csv(input_file[:-4] + "_snp.tsv", sep="\t", index=False)
