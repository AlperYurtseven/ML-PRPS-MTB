import pandas as pd
import sys
import os
from Bio import SeqIO
import numpy as np


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


def add_mutation_info_to_fia(mutation_info_file, fia_file):
    outfile = fia_file + "_MUT_INFO"

    with open(fia_file, "r") as infile:
        fia_lines = infile.readlines()

    with open(mutation_info_file) as snippy_info:
        snippy_info_lines = snippy_info.readlines()

    with open(outfile, "w") as ofile:
        ofile.write(snippy_info_lines[0])
        for line in fia_lines:
            splitted = line.split("\t")
            mutation = splitted[0].strip()
            for line2 in snippy_info_lines:
                splitted2 = line2.split("\t")
                mutation2 = splitted2[0].strip()
                if mutation == mutation2:
                    ofile.write(line2)


def index_finder(bmt, mutation_position):
    # gets binary mutation table and position of mutation, returns strain that has that particular mutation fo feed
    # into snippy
    with open(bmt, 'r') as mutation_file:
        counter = 0
        found_index = None
        mutated_strains = []
        for line in mutation_file.readlines():
            counter += 1
            if counter == 1:
                split_array = line.split('\t')
                for i in range(1, len(split_array) - 1):
                    loc = split_array[i].strip().split(',')[0].strip()[1:]
                    if loc == mutation_position:
                        found_index = i
            else:
                split_array = line.split('\t')
                if split_array[found_index] == '1':
                    mutated_strains.append(split_array[0])

    return mutated_strains[0]


def svm_fia_output_processor(fia_file):
    mutation_positions = []

    with open(fia_file, 'r') as infile:
        lines = infile.readlines()

    cnt = 0
    for line in lines:
        splitted = line.split("\t")
        splitted2 = splitted[1].split(",")
        mutation_positions.append(splitted2[0][1:].strip())
        cnt += 1
        #if cnt == 10:
        #    break

    return mutation_positions


class Mutations:
    """
    Subclass for store the mutation information
    """

    def __init__(self):
        self.pos = ""
        self.type = ""
        self.ref = ""
        self.alt = ""
        self.evidence = ""
        self.ftype = ""
        self.strand = ""
        self.nt_pos = ""
        self.aa_pos = ""
        self.effect = ""
        self.locus_tag = ""
        self.gene = ""
        self.product = ""


def snippy_output_parser(strains, antibiotics):
    stored_mutations = []

    for strain in strains:
        for antibiotic in antibiotics:
            if os.path.isfile(
                    f"./snippy_outputs/{antibiotic}/Resistant/{strain}/snps.tab"):
                file_to_open = f"./snippy_outputs/{antibiotic}/Resistant/{strain}/snps.tab"
                with open(file_to_open, "r") as infile:
                    lines = infile.readlines()

                for line in lines[1:]:
                    temp_mutation = Mutations()
                    splitted = line.split("\t")
                    temp_mutation.pos = splitted[1]
                    temp_mutation.type = splitted[2]
                    temp_mutation.ref = splitted[3]
                    temp_mutation.alt = splitted[4]
                    temp_mutation.evidence = splitted[5]
                    if len(splitted) > 5:
                        temp_mutation.ftype = splitted[6]
                        temp_mutation.strand = splitted[7]
                        temp_mutation.nt_pos = splitted[8]
                        temp_mutation.aa_pos = splitted[9]
                        temp_mutation.effect = splitted[10]
                        temp_mutation.locus_tag = splitted[11]
                        temp_mutation.gene = splitted[12]
                        temp_mutation.product = splitted[13]

                    stored_mutations.append(temp_mutation)

            elif os.path.isfile(
                    f"./snippy_outputs/{antibiotic}/Susceptible/{strain}/snps.tab"):
                file_to_open = f"./snippy_outputs/{antibiotic}/Susceptible/{strain}/snps.tab"
                with open(file_to_open, "r") as infile:
                    lines = infile.readlines()

                for line in lines[1:]:
                    temp_mutation = Mutations()
                    splitted = line.split("\t")
                    temp_mutation.pos = splitted[1]
                    temp_mutation.type = splitted[2]
                    temp_mutation.ref = splitted[3]
                    temp_mutation.alt = splitted[4]
                    temp_mutation.evidence = splitted[5]
                    if len(splitted) > 5:
                        temp_mutation.ftype = splitted[6]
                        temp_mutation.strand = splitted[7]
                        temp_mutation.nt_pos = splitted[8]
                        temp_mutation.aa_pos = splitted[9]
                        temp_mutation.effect = splitted[10]
                        temp_mutation.locus_tag = splitted[11]
                        temp_mutation.gene = splitted[12]
                        temp_mutation.product = splitted[13]

                    stored_mutations.append(temp_mutation)


            elif os.path.isfile(
                   f"./snippy_outputs/{antibiotic}/Resistant/{strain}0/snps.tab"):
                file_to_open = f"./snippy_outputs/{antibiotic}/Resistant/{strain}0/snps.tab"
                with open(file_to_open, "r") as infile:
                    lines = infile.readlines()

                for line in lines[1:]:
                    temp_mutation = Mutations()
                    splitted = line.split("\t")
                    temp_mutation.pos = splitted[1]
                    temp_mutation.type = splitted[2]
                    temp_mutation.ref = splitted[3]
                    temp_mutation.alt = splitted[4]
                    temp_mutation.evidence = splitted[5]
                    if len(splitted) > 5:
                        temp_mutation.ftype = splitted[6]
                        temp_mutation.strand = splitted[7]
                        temp_mutation.nt_pos = splitted[8]
                        temp_mutation.aa_pos = splitted[9]
                        temp_mutation.effect = splitted[10]
                        temp_mutation.locus_tag = splitted[11]
                        temp_mutation.gene = splitted[12]
                        temp_mutation.product = splitted[13]

                    stored_mutations.append(temp_mutation)

            elif os.path.isfile(
                   f"./snippy_outputs/{antibiotic}/Susceptible/{strain}0/snps.tab"):
                file_to_open = f"./snippy_outputs/{antibiotic}/Susceptible/{strain}0/snps.tab"
                with open(file_to_open, "r") as infile:
                    lines = infile.readlines()

                for line in lines[1:]:
                    temp_mutation = Mutations()
                    splitted = line.split("\t")
                    temp_mutation.pos = splitted[1]
                    temp_mutation.type = splitted[2]
                    temp_mutation.ref = splitted[3]
                    temp_mutation.alt = splitted[4]
                    temp_mutation.evidence = splitted[5]
                    if len(splitted) > 5:
                        temp_mutation.ftype = splitted[6]
                        temp_mutation.strand = splitted[7]
                        temp_mutation.nt_pos = splitted[8]
                        temp_mutation.aa_pos = splitted[9]
                        temp_mutation.effect = splitted[10]
                        temp_mutation.locus_tag = splitted[11]
                        temp_mutation.gene = splitted[12]
                        temp_mutation.product = splitted[13]

                    stored_mutations.append(temp_mutation)

            else:
                continue
    return stored_mutations


def all_the_cds_returner(input_file):
    with open(input_file) as infile:
        first_line = infile.readline()
        splitted = first_line.split("\t")
        mutations2 = splitted[1:-6]

    class Protein:
        """
        Subclass for store the protein information
        """

        def __init__(self):
            self.gene = ""
            self.start = 0
            self.end = 0
            self.translation = ""
            self.protein_id = ""
            self.product = ""

    input_gbk = "./reference.gbk"

    recs = [rec for rec in SeqIO.parse(input_gbk, "genbank")]

    all_the_cds = []

    all_the_types = ["CDS", "rRNA", "tRNA"]

    for rec in recs:
        feats = [feat for feat in rec.features if feat.type in all_the_types]

        for feat in feats:
            temp = Protein()
            if 'gene' in feat.qualifiers.keys():
                temp.gene = feat.qualifiers['gene'][0]

            temp.start = feat.location.start
            temp.end = feat.location.end
            try:
                temp.translation = feat.qualifiers['translation'][0]
            except:
                temp.translation = ""
            try:
                temp.protein_id = feat.qualifiers['protein_id'][0]
            except:
                temp.protein_id = ""
            try:
                temp.product = feat.qualifiers['product'][0]
            except:
                temp.product = ""
            try:
                temp.gene = feat.qualifiers['gene'][0]
            except:
                temp.gene = ""
            try:
                temp.note = feat.qualifiers['note'][0]
            except:
                temp.note = ""
            all_the_cds.append(temp)

    return all_the_cds, mutations2


def svm_fia_mutation_information_adder(fia_file, mutations, dataframe):
    card_db = {
        "katG": "Prothionamide & Isoniazid resistance\tFirst & Second Line\t- & Nicotinamide derivative",
        "rpoC": "Rifampicin resistance\tFirst line\t-",
        "fusA2": "Fusidic acid resistance\t-\t-",
        "rpoB": "Rifampicin resistance\tFirst line\t-",
        "ubiA": "Ethambutol resistance\tFirst line\t-",
        "ponA1": "Rifampicin resistance\tFirst line\t-",
        "mmaA4": "Isoniazid resistance\tFirst line\tNicotinamide derivative",
        "embB": "Ethambutol & Rifampicin resistance\tFirst line\t- & -",
        "ndhA": "Isoniazid resistance\tSecond line\tNicotinamide derivative",
        "rpsA": "Pyrazinamide resistance\tFirst line\tNicotinamide derivative",
        "ahpD": "Isoniazid resistance\tFirst line\tNicotinamide derivative",
        "eccD2": "eccC5 & eccB5 = Fluoroquinolone resistance\tSecond line\tFluoroquinolones",
        "gyrA": "Fluoroquinolone resistance\tSecond line\tFluoroquinolones",
        "mshC": "Ethionamide & Isoniazid resistance\tFirst & Second Line\tNicotinamide derivative & -",
        "rpsL": "Streptomycin resistance\tFirst line\tAminoglycosides",
        "ileS": "Mupirocin resistance\t-\t-",
        "aftA": "Ethambutol resistance\tFirst line\t-",
        "hsaA": "Halobacterium salinarum Pactamycin resistance\t-\t-",
        "sigE": "sigI Isoniazid resistance\tFirst line\tNicotinamide derivative",
        "eccA2": "eccC = Fluoroquinolone resistance\tSecond line\tFluoroquinolones",
        "eccE3": "eccC = Fluoroquinolone resistance\tSecond line\tFluoroquinolones",
        "embA": "Ethambutol & Rifampicin resistance\tFirst line\t-",
        "proA": "Isoniazid resistance\tFirst line\tNicotinamide derivative",
        "fabG1": "Ethionamide & Isoniazid resistance\tFirst & Second Line\tNicotinamide derivative",
        "inhA": "Isoniazid resistance\tFirst line\tNicotinamide derivative",
        "rrs": "Aminoglycoside resistance\tFirst & Second Line\tAminoglycosides",
        "eis": "Kanamycin resistance\tSecond line\tAminoglycosides",
        "gyrB": "Fluoroquinolone resistance\tSecond line\tFluoroquinolones",
        "ethA": "Ethionamide resistance\tSecond line\tNicotinamide derivative",
        "ethR": "Ethionamide resistance\tSecond line\tNicotinamide derivative",
        "tlyA": "Aminoglycoside resistance\tFirst & Second Line\tAminoglycosides",
        "gidB": "Aminoglycoside resistance\tFirst & Second Line\tAminoglycosides",
        "ahpC": "Isoniazid resistance\tFirst line\tNicotinamide derivative",
        "ndh": "Isoniazid resistance\tFirst line\tNicotinamide derivative",
        "embC": "Ethambutol resistance\tFirst line\t-",
        "embR": "Ethambutol resistance\tFirst line\t-",
        "iniA": "Ethambutol resistance\tFirst line\t-",
        "iniC": "Ethambutol resistance\tFirst line\t-",
        "manB": "Ethambutol resistance\tFirst line\t-",
        "rmlD": "Ethambutol resistance\tFirst line\t-"

    }

    mutation_positions = []

    all_the_cds, mutations2 = all_the_cds_returner(dataframe)

    with open(fia_file, 'r') as infile:
        lines = infile.readlines()

    out_lines = []
    cnt = 0
    out_lines.append(
        "CHROM\t(POS,REF,ALT,MTYPE)\tTYPE\tNT_POS\tAA_POS\tEFFECT\tGENE\tPRODUCT\tCARD_DB\tLINE_OF_DRUG\tPHARMA_GROUP\tDISTANCE")
    for line in lines:
        mut_appended = False
        splitted = line.split("\t")
        splitted2 = splitted[1].split(",")
        pos = splitted2[0][1:].strip()
        mutation_positions.append(splitted2[0][1:].strip())
        for mut in mutations:
            if mut.pos == pos:
                if mut.gene in card_db.keys():
                    out_lines.append(
                        line.strip() + "\t" + mut.ftype + "\t" + mut.nt_pos + "\t" + mut.aa_pos + "\t" + mut.effect + "\t" + mut.gene + "\t" + mut.product.strip() + "\t" +
                        card_db[mut.gene])
                    mut_appended = True
                else:
                    out_lines.append(
                        line.strip() + "\t" + mut.ftype + "\t" + mut.nt_pos + "\t" + mut.aa_pos + "\t" + mut.effect + "\t" + mut.gene + "\t" + mut.product.strip())
                    mut_appended = True
                break

        if not mut_appended:
            closest_j = None
            closest_val = 9999999
            for j in all_the_cds:
                if abs(j.start - int(pos)) < closest_val:
                    closest_j = j
                    closest_val = abs(j.start - int(pos))

                if abs(j.end - int(pos)) < closest_val:
                    closest_j = j
                    closest_val = abs(j.end - int(pos))

            if closest_j.gene in card_db.keys():
                out_lines.append(line.strip() + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + str(
                    closest_j.gene) + "\t" + str(closest_j.product) + "\t" + card_db[
                                     closest_j.gene] + "\t" + "Closest_dist: " + str(
                    closest_val))
            else:
                out_lines.append(line.strip() + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + "" + "\t" + str(
                    closest_j.gene) + "\t" + str(closest_j.product) + "\t\t\t\t" + "Closest_dist: " + str(
                    closest_val))

        cnt += 1

    with open(fia_file + "_mutINFO.tsv", "w") as ofile:
        for line in out_lines:
            ofile.write(str(line) + "\n")


def mutation_file_creator(dataframe, outfile):

    with open(outfile, "w") as ofile:
        with open(dataframe) as infile:
            line = infile.readline
            splitted = line.split("\t")
            for mut in splitted:
                if mut.startswith("("):
                    ofile.write(f"1\t{mut}\n")


def annotation_file_creator(dataframe, mutation_list):

    mutation_file_creator(dataframe, mutation_list)

    mut_positions = svm_fia_output_processor(mutation_list)

    mutated_strains = []
    for mut_pos in mut_positions:
        mutated_strains.append(index_finder(
            dataframe,
            mut_pos))

    mutations = snippy_output_parser(mutated_strains)
    svm_fia_mutation_information_adder(
        mutation_list, mutations, dataframe)


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


def snippy_runner_reference(main_path, list_of_strains, reference, output_dir):
    set_of_strains = set(list_of_strains)

    main_path = main_path

    for strain in set_of_strains:
        input_path = main_path + strain
        os.system(f"snippy --outdir {output_dir}/{strain[:-4]}_snippy --ref {reference} --ctgs {input_path}")

def snippy_runner(main_path, list_of_strains, output_dir):
    set_of_strains = set(list_of_strains)

    main_path = main_path

    for strain in set_of_strains:
        input_path = main_path + strain
        os.system(f"snippy --outdir {output_dir}/{strain[:-4]}_snippy --ctgs {input_path}")
