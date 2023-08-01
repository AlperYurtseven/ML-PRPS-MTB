import pandas as pd
import sklearn
from sklearn.metrics import classification_report
from sklearn.svm import SVC


def full_pipeline(percentage):
    
    print("%s pipeline started" % percentage)

    print("%s percentage column extractor" % percentage)

    #percentage_column_extractor(percentage, "./20032023_sonia_results.tsv")

    print("%s columns dropper" % percentage)

    #columns_dropper(percentage)

    columns_dropper_order(percentage)

    print("%s seperate antibiotic creator" % percentage)

    #seperate_antibiotic_creator("/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_%s.tsv" % percentage)

    #seperate_antibiotic_creator("/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_%s.tsv" % percentage)

    seperate_antibiotic_creator(
        "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_ordered_%s.tsv" % percentage)

    antibiotics = ["amikacin", "capreomycin", "ethionamide", "kanamycin", "ofloxacin", "streptomycin"]
    for abiotic in antibiotics:

        print("%s svm model training for: %s" % (percentage, abiotic))

        # svm_ml(
        #     "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_%s_%s.tsv" % (percentage,abiotic))

        svm_ml(
            "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_ordered_%s_%s.tsv" % (
            percentage, abiotic))

        print("%s snippy info adding for: %s" % (percentage, abiotic))

        # add_snippy_info_to_fia(
        # "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_%s_%s_SVM_FIA_TOP10" % (percentage,abiotic))

        add_snippy_info_to_fia(
            "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_ordered_%s_%s_SVM_FIA_TOP20" % (
            percentage, abiotic))


def max_value_calculator(input_file):

    max_score = 0.0

    with open(input_file, "r") as infile:
        lines = infile.readlines()

    for line in lines:
        splitted = line.split("\t")

        if len(splitted) < 1:
            continue
        if splitted[1].strip() == "score":
            continue
        if splitted[1].strip() == "":
            continue
        else:
            phy_score = splitted[1].strip()
            phy_score= phy_score.replace(",", ".")
            #print(splitted[0])
            phy_score = float(phy_score)
            if phy_score > max_score:
                max_score = phy_score
    return max_score


def percentage_column_extractor(percentage, input_file):

    max_value = max_value_calculator(input_file)

    percentage = int(percentage)

    threshold = (max_value*percentage) / 100

    #print(threshold)

    out_file = input_file[:-4] + "_less_than_threshold_%s" % str(percentage)
    out_file2 = input_file[:-4] + "_more_than_threshold_%s" % str(percentage)

    with open(input_file, "r") as infile:
        lines = infile.readlines()

    with open(out_file, "w") as ofile:
        with open(out_file2, "w") as ofile2:
            for line in lines:
                splitted = line.split("\t")
                if len(splitted) < 1:
                    continue
                if splitted[1].strip() == "score":
                    continue
                if splitted[1].strip() == "":
                    continue
                else:
                    phy_score = splitted[1].strip()
                    phy_score= phy_score.replace(",", ".")
                    phy_score = float(phy_score)
                    if phy_score < threshold:
                        ofile.write(line.strip())
                        ofile.write("\n")
                    if phy_score >= threshold:
                        ofile2.write(line.strip())
                        ofile2.write("\n")


def columns_dropper(percentage):

    input_file = "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected.tsv"

    df = pd.read_csv(input_file, sep="\t", dtype={"name/position": str, "outcome": int, "*": int})

    mutations_to_drop = []
    with open("./20032023_sonia_results_less_than_threshold_%s" % percentage) as infile:
        lines = infile.readlines()

        for line in lines:
            splitted = line.split("\t")
            mutations_to_drop.append(splitted[0].strip())


    dropped_df = df.drop(mutations_to_drop, axis=1)

    dropped_df.to_csv(input_file[:-4] + "_after_pyhlogeny_%s.tsv" % percentage, sep="\t", index=False)


def add_snippy_info_to_fia(fia_file, annotation_file):

    outfile = fia_file + "_MUT_INFO"

    with open(fia_file, "r") as infile:
        fia_lines = infile.readlines()

    #with open("./20032023_sonia_results_MUT_INFO.tsv") as snippy_info:
    with open(annotation_file) as snippy_info:
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


def top_percent_dropper(percent, annotation_file_ordered):
    top_mutations = []
    mutations_to_drop = []

    #with open("./20032023_sonia_results_MUT_INFO_ORDERED.tsv") as infile:
    with open(annotation_file_ordered) as infile:
        lines = infile.readlines()

    line_count = len(lines)-1

    percent = int(percent)

    limit = percent * line_count / 100
    limit = int(limit)

    i = 0

    for line in lines:
        if i >= limit:
            splitted = line.split("\t")
            mutations_to_drop.append(splitted[0])
            i += 1
        else:
            i += 1

    print(str(len(mutations_to_drop)) + " / " + str(len(lines)-1))

    return mutations_to_drop


def columns_dropper_order(percentage, input_file):

    #input_file = "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_column_corrected.tsv"

    df = pd.read_csv(input_file, sep="\t", dtype={"name/position": str, "outcome": int, "*": int})

    mutations_to_drop = top_percent_dropper(percentage)

    dropped_df = df.drop(mutations_to_drop, axis=1)

    dropped_df.to_csv(input_file[:-4] + "_after_pyhlogeny_ordered_%s.tsv" % percentage, sep="\t", index=False)



if __name__ == "__main__":

    percentages = ["30"]

    for perc in percentages:
        full_pipeline(perc)