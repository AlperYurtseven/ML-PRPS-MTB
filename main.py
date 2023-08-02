import sys
import os
from functions import snippy_runner_reference, binary_table_combiner, filter_percentage, seperate_antibiotic_creator, zero_col_dropper, annotation_file_creator, drop_non_snp_columns, columns_dropper_order
from binary_mutation_table_creator_manual import RawDataProcessor
from ml_svm import svm_ml
from ml_rf import rf_auto_ml
from add_mutation_info_to_fia import add_mutation_info_to_fia
from ml_svm_random_drop import svm_ml_random_drop
from ml_rf_random_drop import rf_auto_ml_random_drop


def main(data_path, reference_path, snippy_out_path, binary_table_path, filter_percentage_value, phylogeny_percentage, file_extension=".fna"):

    # Run Snippy for variant calling

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

    # Process variant calling output
     
    for antibiotic in antibiotics_list:
        RawDataProcessor(f"{snippy_out_path}/{antibiotic}/Resistant", f"{snippy_out_path}/{antibiotic}/Susceptible", f"{binary_table_path}/{antibiotic}_binary_table.tsv")

    binary_table_combiner(f"{binary_table_path}", f"combined_binary_table", antibiotics_list)
    
    filter_percentage(f"{binary_table_path}/combined_binary_table.tsv", filter_percentage_value, antibiotics_list, f"{binary_table_path}/combined_binary_table")

    drop_non_snp_columns(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}.tsv")

    seperate_antibiotic_creator(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_snp.tsv")

    for antibiotic in antibiotics_list:
        zero_col_dropper(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_snp_{antibiotic}.tsv")

    annotation_file_creator(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}.tsv", f"{binary_table_path}/mutation_annotation_file")
    
    # Model training before PRPS score training 

    for antibiotic in antibiotics_list:
        svm_ml(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}_snp_dropped_zero_cols.tsv")
        rf_auto_ml(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}_snp_dropped_zero_cols.tsv")
    
    #TODO PRPS score calculation

    # Model training after PRPS score training same procedures/filters are applied

    columns_dropper_order(phylogeny_percentage, f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_snp.tsv")

    seperate_antibiotic_creator(
        f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_snp_after_pyhlogeny_ordered_{phylogeny_percentage}.tsv")
    
    for antibiotic in antibiotics_list:
        zero_col_dropper(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_snp_after_pyhlogeny_ordered_{phylogeny_percentage}_{antibiotic}.tsv")
    
    for antibiotic in antibiotics_list:
        svm_ml(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_snp_after_pyhlogeny_ordered_{phylogeny_percentage}_{antibiotic}_dropped_zero_cols.tsv")
        rf_auto_ml(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_snp_after_pyhlogeny_ordered_{phylogeny_percentage}_{antibiotic}_dropped_zero_cols.tsv")


    # Train same amount of randomly deleted feature models part to check PRPS working as expected
    
    for antibiotic in antibiotics_list:
        svm_ml_random_drop(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}_snp_dropped_zero_cols.tsv", phylogeny_percentage)
        rf_auto_ml_random_drop(f"{binary_table_path}/combined_binary_table_{filter_percentage_value}_{antibiotic}_snp_dropped_zero_cols.tsv", phylogeny_percentage)
    

if __name__ == "__main__":
    main("./data/", "./reference.gbff", "./snippy_outputs", "./binary_tables", 0.2, 30)