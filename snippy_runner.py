import os


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

# all_strains = os.listdir("/scratch/SCRATCH_SAS/PATRIC_AMR_data/Mycobacterium/zcombinations/")

# strains_to_run_snippy = []
# for i in all_strains:
#     if i.endswith(".fna"):
#         strains_to_run_snippy.append(i)

# snippy_runner(strains_to_run_snippy)
