import pandas as pd
import numpy as np
import sys
import os


def binary_table_combiner(binary_table_path, output_name, antibiotics):

    files_list = os.listdir(binary_table_path)
    files = []
    for file in files_list:
        files.append(f"{binary_table_path}/{file}")

    all_the_mutations = []
    all_the_strains = []

    frames = []

    count = 0
    for f in files:
        df = pd.read_csv(f, sep="\t", dtype={"name/position": str, "outcome": str, "*": str}, engine="pyarrow")
        df = df.astype(str)
        for i in df.columns:
            if i == "name/position":
                continue
            elif i == "outcome":
                continue
            else:
                if i not in all_the_mutations:
                    all_the_mutations.append(i)

    all_df_dicts = []
    for f in files:
        df = pd.read_csv(f, sep="\t", dtype={"name/position": str, "outcome": str, "*": str})
        df = df.astype(str)
        df_dict = df.to_dict(orient="index")

        for antibiotic in antibiotics:
            if antibiotic in f:
                for i in df_dict.keys():
                    df_dict[i][antibiotic] = df_dict[i]["outcome"]
                    df_dict[i].pop("outcome")

        temp_dict = {}

        for key in df_dict.keys():
            temp_name = df_dict[key]["name/position"]
            df_dict[key].pop("name/position")
            temp_dict[temp_name] = df_dict[key]
            all_the_strains.append(temp_name)

        all_df_dicts.append(temp_dict)
    # print("a")

    combination_dict = {}

    all_the_strains_set = set(all_the_strains)

    for strain in all_the_strains_set:
        temp_dict = {}
        for dic in all_df_dicts:
            if strain in dic.keys():
                for key2 in dic[strain]:
                    temp_dict[key2] = dic[strain][key2]

        for mut in all_the_mutations:
            if mut in temp_dict.keys():
                continue
            else:
                temp_dict[mut] = 0
        combination_dict[strain] = temp_dict

    df2 = pd.DataFrame.from_dict(combination_dict, orient="index")

    df2 = df2.astype(str)
    cols = df2.columns.tolist()

    for antibiotic in antibiotics:
        cols.remove(antibiotic)

    for antibiotic in antibiotics:
        cols.append(antibiotic)

    df2 = df2[cols]

    # print(df2)
    # df2 = df2.replace(pd.NA, 0)
    # print(df2)
    # df2 = df2.replace("nan", value="1")
    df2 = df2.astype(str)
    # df2.to_csv("./deneme_out_combination.tsv", sep="\t")
    df2.to_csv(f"{binary_table_path}/{output_name}.tsv", sep="\t")

