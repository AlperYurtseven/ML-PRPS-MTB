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

