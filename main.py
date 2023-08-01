import sys
import os
from functions import snippy_runner_reference, binary_table_combiner, filter_percentage, seperate_antibiotic_creator, zero_col_dropper, annotation_file_creator, drop_non_snp_columns
from binary_mutation_table_creator_manual import RawDataProcessor
from ml_svm import svm_ml
from ml_rf import rf_auto_ml
from add_mutation_info_to_fia import add_mutation_info_to_fia


def main(data_path, reference_path, snippy_out_path, binary_table_path, filter_percentage_value, file_extension=".fna"):

    antibiotics_list = os.listdir(data_path)

    for antibiotic in antibiotics_list:
        resistance_status = os.listdir(f"{data_path}/{antibiotic}")

        for rs in resistance_status:
            temp_strain_list = os.listdir(f"{data_path}/{antibiotic}/{rs}")
            strains_to_run_snippy = []
            for s in temp_strain_list:
                if s.endswith(file_extension):
                    strains_to_run_snippy.append(s)

            snippy_runner_reference(f"{data_path}/{antibiotic}/{rs}/", strains_to_run_snippy, reference_path, f"{snippy_out_path}/{antibiotic}/{rs}")


    for antibiotic in antibiotics_list:
        RawDataProcessor(f"{snippy_out_path}/{antibiotic}/Resistant", f"{snippy_out_path}/{antibiotic}/Susceptible", f"{binary_table_path}/{antibiotic}_binary_table.tsv")

    binary_table_combiner(f"{binary_table_path}", f"combined_binary_table", antibiotics_list)
    
    filter_percentage(f"{binary_table_path}/combined_binary_table.tsv", filter_percentage_value, antibiotics_list, f"{binary_table_path}/combined_binary_table")

    seperate_antibiotic_creator(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}.tsv")

    for antibiotic in antibiotics_list:
        drop_non_snp_columns(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}.tsv")
        zero_col_dropper(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}_snp.tsv")

    annotation_file_creator(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}.tsv", f"{binary_table_path}/mutation_annotation_file")
    
    for antibiotic in antibiotics_list:
        svm_ml(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}_snp_dropped_zero_cols.tsv.tsv")
        rf_auto_ml(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}_snp_dropped_zero_cols.tsv.tsv")
        



if __name__ == "__main__":
    main("./data/", "./reference.gbff", "./snippy_outputs", "./binary_tables", 0.2)