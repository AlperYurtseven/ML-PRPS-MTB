import sys
import pandas as pd


def bmf_to_vcf(input_df, opath):
    pos_dict = {}
    ref_allele_dict = {}
    alt_allele_dict = {}
    type_dict = {}
    strains_dict = {}
    all_the_strains = []
    with open(input_df) as infile:
        lines = infile.readlines()
    for line in lines:
        splitted = line.split("\t")
        if splitted[0] == "name/position":
            continue
        else:
            all_the_strains.append(splitted[0])

    cnt = 0
    for i in all_the_strains:
        strains_dict[cnt] = i
        cnt += 1

    strains_dict_other = {}
    for i in strains_dict.keys():
        strains_dict_other[strains_dict[i]] = i

    df = pd.read_csv(input_df, sep="\t", dtype={"name/position": str, "outcome": int, "*": int}, low_memory=False)

    #print(df.shape)
    c = df.iloc[0][1]
    d = df.shape
    flag = True
    with open(opath, "w") as ofile:
        ofile.write("##fileformat=VCFv4.2\n")
        ofile.write("##FILTER=<ID=PASS,Description=""All filters passed"">\n")
        ofile.write("##fileDate=20211216\n")
        ofile.write("##source=freeBayes v1.3.2-dirty\n")
        ofile.write("##reference=reference/ref.fa\n")
        ofile.write("##contig=<ID=1,length=4411532>\n")
        ofile.write("##phasing=none\n")
        ofile.write("##commandline=""freebayes -p 2 -P 0 -C 2 -F 0.05 --min-coverage 10 --min-repeat-entropy 1.0 -q 13 -m 60 --strict-vcf -f reference/ref.fa snps.bam --region NC_000962.3:0-299008""\n")
        ofile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=""Total read depth at the locus"">\n")
        ofile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">\n")
        ofile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for i in strains_dict.keys():
            ofile.write("\t" + strains_dict[i])
        ofile.write("\n")
        with open(input_df) as infile:
            lines = infile.readlines()
            for line in lines:
                splitted = line.split("\t")
                index = 0
                if line.startswith("name"):
                    for s in splitted:
                        if s.strip() == "name/position":
                            continue
                        if s.strip() == "outcome":
                            continue
                        else:
                            # a = s.split(",")
                            pos = int(s.split(",")[0][1:])
                            ref_all = s.split(",")[1][2]
                            alt_all = s.split(",")[2][0]
                            type = s.split(",")[3][2:-2]

                            pos_dict[index] = pos
                            ref_allele_dict[index] = ref_all
                            alt_allele_dict[index] = alt_all
                            type_dict[index] = type
                            index += 1
                elif flag:
                    flag = False
                    threshold = df.shape[1] - 2
                    for i in range(1,threshold):
                        ofile.write(
                            "1" + "\t" + str(pos_dict[i]) + "\t" + "." + "\t" + str(ref_allele_dict[i]) + "\t" + str(
                                alt_allele_dict[i]) + "\t" + "." + "\t" + "." + "\t" + "PR" + "\t" + "GT")
                        for j in strains_dict.keys():
                            c = df.iloc[j][i]
                            if str(df.iloc[j][i]) == "1":
                                ofile.write("\t" + "1/1")
                            elif str(df.iloc[j][i]) == "0":
                                ofile.write("\t" + "0/0")
                            else:
                                ofile.write("\t" + "./.")
                        ofile.write("\n")


#antibiotics = ["amikacin", "capreomycin", "ethionamide", "kanamycin", "ofloxacin", "streptomycin"]

bmf_to_vcf(sys.argv[1], sys.argv[2])

# for abiotic in antibiotics:
#     bmf_to_vcf("/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_%s_name_corrected.tsv" % abiotic, "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_%s_vcf.vcf" % abiotic)