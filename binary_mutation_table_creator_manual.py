import sys
import pandas as pd
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from datetime import datetime
import argparse

class RawDataProcessor:
    """
    Class for process raw data takes resistance and susceptible data paths as well as binary output path
    Creates mutation binary file from snippy output.

    """

    def __init__(self, path_for_resistance_data, path_for_susceptible_data, binary_output_path):

        self.main(path_for_resistance_data, path_for_susceptible_data, binary_output_path)

    class snps_class():
        """
       Class to store the snps with their chromosome, position, reference allele, alternative allele,
       quality and mutation type
       """
        chrom = ""
        pos = ""
        ref = ""
        alt = ""
        qual = 0.00
        mut_type = ""

    def main(self, path_for_resistance_data, path_for_susceptible_data, binary_output_path):

        """
        Main Function to run and create binary mutation file

        Parameters
        ----------
        path_for_resistance_data : str
            path to resistance data of specific antibiotic that has snippy outputs stored in the folders.

        path_for_susceptible_data : str
            path to susceptible data of specific antibiotic that has snippy outputs stored in the folders.

        binary_output_path : str
            path to binary output, will be created during runtime, if it is already exist, will be overwritten

        Return
        ----------

        None

        """
        path = path_for_resistance_data
        path2 = path_for_susceptible_data
        output_path = binary_output_path

        all_directories_res = []
        all_directories_sus = []

        for i in os.listdir(path):
            if i.endswith("_snippy"):
                all_directories_res.append(i)

        for i in os.listdir(path2):
            if i.endswith("_snippy"):
                all_directories_sus.append(i)

        all_the_mutations = {}

        file_name_and_pos_dict_res = {}
        file_name_and_pos_dict_sus = {}

        all_the_tuples = []

        """
        Part to remove directories from list if snps.vcf file did not created for them by snippy

        Ensures that snps.vcf will be exist on the paths will be open in next part
        """

        for b in range(len(all_directories_res)):
            try:
                new_dir = path + all_directories_res[b] + "/snps.vcf"
                with open(new_dir) as infile:
                    lines = infile.readlines()

                snp_class_list = []

                for line in lines:
                    if line[0] == '#':
                        continue
                    else:
                        splitted = line.split("\t")
                        temp = self.snps_class()
                        temp.chrom = splitted[0]
                        temp.pos = splitted[1]
                        temp.ref = splitted[3]
                        temp.alt = splitted[4]
                        temp.qual = splitted[5]
                        splitted_info = splitted[7].split(';')
                        for p in splitted_info:
                            if p.startswith("TYPE"):
                                splitted_info2 = p.split('=')
                                temp.mut_type = splitted_info2[1]

                        snp_class_list.append(temp)

                for a in snp_class_list:
                    temp_str = a.ref + "," + a.alt
                    temp_tuple = (a.pos, temp_str, a.mut_type)
                    # all_the_tuples2.append(temp_tuple)
                    if temp_tuple not in all_the_tuples:
                        all_the_tuples.append(temp_tuple)

            except FileNotFoundError:
                continue



        for c in range(len(all_directories_sus)):
            try:
                new_dir = path2 + all_directories_sus[c] + "/snps.vcf"
                with open(new_dir) as infile:
                    lines = infile.readlines()

                snp_class_list = []

                for line in lines:
                    if line[0] == '#':
                        continue
                    else:
                        splitted = line.split("\t")
                        temp = self.snps_class()
                        temp.chrom = splitted[0]
                        temp.pos = splitted[1]
                        temp.ref = splitted[3]
                        temp.alt = splitted[4]
                        temp.qual = splitted[5]
                        splitted_info = splitted[7].split(';')
                        for p in splitted_info:
                            if p.startswith("TYPE"):
                                splitted_info2 = p.split('=')
                                temp.mut_type = splitted_info2[1]

                        snp_class_list.append(temp)

                for a in snp_class_list:
                    temp_str = a.ref + "," + a.alt
                    temp_tuple = (a.pos, temp_str, a.mut_type)
                    # all_the_tuples2.append(temp_tuple)
                    if temp_tuple not in all_the_tuples:
                        all_the_tuples.append(temp_tuple)
            except FileNotFoundError:
                continue

        all_the_tuples_int = []

        for tup in all_the_tuples:
            new_tup = (int(tup[0]), tup[1], tup[2])
            all_the_tuples_int.append(new_tup)

        all_the_tuples_int.sort(key=lambda x: x[0])

        cnt = 0
        for tup in all_the_tuples_int:
            all_the_mutations[cnt] = tup
            cnt += 1

        for c in range(len(all_directories_sus)):
            try:
                new_dir = path2 + all_directories_sus[c] + "/snps.vcf"

                with open(new_dir) as infile:
                    lines = infile.readlines()

                snp_class_list = []

                for line in lines:
                    if line[0] == '#':
                        continue
                    else:
                        splitted = line.split("\t")
                        temp = self.snps_class()
                        temp.chrom = splitted[0]
                        temp.pos = splitted[1]
                        temp.ref = splitted[3]
                        temp.alt = splitted[4]
                        temp.qual = splitted[5]
                        splitted_info = splitted[7].split(';')
                        for p in splitted_info:
                            if p.startswith("TYPE"):
                                splitted_info2 = p.split('=')
                                temp.mut_type = splitted_info2[1]

                        snp_class_list.append(temp)

                temp_list = [0] * len(all_the_mutations)

                for a in snp_class_list:
                    temp_str = a.ref + "," + a.alt
                    temp_tuple = (int(a.pos), temp_str, a.mut_type)

                    if temp_tuple in all_the_mutations.values():
                        x = list(all_the_mutations.keys())[list(all_the_mutations.values()).index(temp_tuple)]
                        temp_list[x] = 1

                file_name_temp = all_directories_sus[c][:-11]

                file_name_and_pos_dict_sus[file_name_temp] = temp_list

            except FileNotFoundError:
                continue

        for c in range(len(all_directories_res)):
            try:
                new_dir = path + all_directories_res[c] + "/snps.vcf"

                with open(new_dir) as infile:
                    lines = infile.readlines()

                snp_class_list = []

                for line in lines:
                    if line[0] == '#':
                        continue
                    else:
                        splitted = line.split("\t")
                        temp = self.snps_class()
                        temp.chrom = splitted[0]
                        temp.pos = splitted[1]
                        temp.ref = splitted[3]
                        temp.alt = splitted[4]
                        temp.qual = splitted[5]
                        splitted_info = splitted[7].split(';')
                        for p in splitted_info:
                            if p.startswith("TYPE"):
                                splitted_info2 = p.split('=')
                                temp.mut_type = splitted_info2[1]

                        snp_class_list.append(temp)

                temp_list = [0] * len(all_the_mutations)

                for a in snp_class_list:
                    temp_str = a.ref + "," + a.alt
                    temp_tuple = (int(a.pos), temp_str, a.mut_type)

                    if temp_tuple in all_the_mutations.values():
                        x = list(all_the_mutations.keys())[list(all_the_mutations.values()).index(temp_tuple)]
                        temp_list[x] = 1

                file_name_temp = all_directories_res[c][:-11]

                file_name_and_pos_dict_res[file_name_temp] = temp_list

            except FileNotFoundError:
                continue

        with open(output_path, "w") as out_file:

            out_file.write("name/position")
            for t in range(len(all_the_mutations)):
                out_file.write("\t")
                out_file.write(str(all_the_mutations[t]))

            out_file.write("\t")
            out_file.write("outcome")
            out_file.write("\n")

            for k in file_name_and_pos_dict_res.keys():
                out_file.write(str(k))
                for e in file_name_and_pos_dict_res[k]:
                    out_file.write("\t")
                    out_file.write(str(e))

                out_file.write("\t")
                out_file.write(str(1))
                out_file.write("\n")

            for k in file_name_and_pos_dict_sus.keys():
                out_file.write(str(k))
                for e in file_name_and_pos_dict_sus[k]:
                    out_file.write("\t")
                    out_file.write(str(e))

                out_file.write("\t")
                out_file.write(str(0))
                out_file.write("\n")

        out_file.close()


def main():

    RawDataProcessor(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()