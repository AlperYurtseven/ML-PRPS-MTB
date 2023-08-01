import pandas as pd


def filter_percentage(file, percentage, antibiotics, output_file):

    dataframe = pd.read_csv(file, sep="\t", dtype={"name/position": str, "outcome": str, "*": str}, low_memory=False)

    total_rows = dataframe.shape[0]

    threshold = (total_rows / 100) * percentage

    threshold = int(threshold)

    cols_to_drop = []

    for col in dataframe.columns:
        if col in antibiotics:
            continue
        elif col == "name/position":
            continue
        elif col == "outcome":
            continue
        elif dataframe[col].sum() < threshold:
            cols_to_drop.append(col)

    
    return_df = dataframe.drop(cols_to_drop, axis=1)
    
    return_df.to_csv(f"{output_file}_{str(percentage)}.tsv", sep="\t")

    #print("Amount of cols dropped with %.1f percent = %d out of %d" % (percentage, len(cols_to_drop), dataframe.shape[1]-7))
