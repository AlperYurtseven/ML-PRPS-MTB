import pandas as pd


input_file = "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected.tsv"

df = pd.read_csv(input_file, sep="\t", dtype={"name/position": str, "outcome": int, "*": int})

mutations_to_drop = []
with open("./sonia_results_less_than_threshold_35") as infile:
    lines = infile.readlines()

    for line in lines:
        splitted = line.split("\t")
        mutations_to_drop.append(splitted[0].strip())


dropped_df = df.drop(mutations_to_drop, axis=1)

dropped_df.to_csv(input_file[:-4] + "_after_pyhlogeny_35.tsv", sep="\t", index=False)