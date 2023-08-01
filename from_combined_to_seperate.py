import pandas as pd

def seperate_antibiotic_creator(file, antibiotics):
    df = pd.read_csv(file, sep="\t", dtype={"name/position": str, "outcome": str, "*": str}, low_memory=False, index_col=0)

    #df["amikacin"] = df["amikacin"].astype("int")

    for antibiotic in antibiotics:
        temp_list = [item for item in antibiotics if item not in antibiotic]
        df[antibiotic] = df[antibiotic].replace(0.0, "0")
        df[antibiotic] = df[antibiotic].replace(1.0, "1")

        df = df.fillna(value="?")

        df.to_csv(file[:-4] + "_filled.tsv", sep="\t", index=False)

        antibiotic_df = df.drop(temp_list, axis=1)
        antibiotic_df = antibiotic_df[antibiotic_df.antibiotic != "?"]
        new_name = f"{file[:-4]}_{antibiotic}.tsv"
        antibiotic_df.rename(columns={antibiotic:"outcome"}, inplace=True)
        antibiotic_df.to_csv(new_name, sep="\t", index=False)
