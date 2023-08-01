import os
import sys
import pandas as pd
from Bio import SeqIO


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

   