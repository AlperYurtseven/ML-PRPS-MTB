import pandas as pd
import os
import seaborn as sns



def load_renaming_and_features(bact_dataloc, datafile, featurefile, origformat=False):
    """ Matches strain identifiers of Panacota output, 

    Parameters
    ----------
    bact_dataloc : str, this is where the data is located
    datafile : str, name of file panacota produce with the renames strain names
    featurefile : merged feature full file name path, result of load_features_to_table func
    origformat (bool): A flag that indicates whether the data file is in its original format. Default is False.
    
    Returns:

    id_match (pd.DataFrame): A DataFrame that contains the matched strain identifiers.


    """
    ecoli_strains=pd.read_csv(featurefile, sep="\t", index_col="Unnamed: 0")
    
    if "name/position" in ecoli_strains.columns:
        ecoli_strains = ecoli_strains.rename(columns={"name/position":"strain"})
    
    ecoli_strains.loc[:,"strain"]=ecoli_strains.strain.astype("str")
    
    if origformat == False:
        ecoli_renamed=pd.read_csv(f"{bact_dataloc}/{datafile}",
                                  sep="\t", header=None)
        orig_strain=ecoli_renamed[1].str.split("/", expand=True)[2].str.split(".fna_p", expand=True)[0]
        ecoli_renamed[1]=orig_strain
        ecoli_renamed=ecoli_renamed[[0,1]]
        ecoli_renamed=ecoli_renamed.rename(columns={0:"panacota_renamed",1:"strain"})
        
    else:
        
        ecoli_renamed=pd.read_csv(f"{bact_dataloc}/{datafile}",
                                  sep="\t")
        orig_strain=ecoli_renamed.iloc[:,1].str.split("/", expand=True)[3].str.split(".fna_p", expand=True)[0]
        ecoli_renamed.iloc[:,1]=orig_strain
        ecoli_renamed=ecoli_renamed.iloc[:,[0,1]]
        ecoli_renamed=ecoli_renamed.rename(columns={"gembase_name":"panacota_renamed","orig_name":"strain"})

    id_match=ecoli_renamed.merge(ecoli_strains, on="strain")
    return id_match

def strain_colors(numcolors=14,col_palette="viridis"):
    #sets pallette
    #numcolors = len(id_match.drug.unique())
    #colors = list(reversed(sns.color_palette("viridis", numcolors).as_hex()))
    colors = sns.color_palette(col_palette, numcolors).as_hex()
    return colors

def write_color_to_file(
    id_match,
    antibiotic,
    color,
    bact_name,
    melted=True):
    
    """
    Parameters:

    id_match (pd.DataFrame): A DataFrame that contains the matched strain identifiers.
    antibiotic (str): The name of the antibiotic.
    color (str): The hex code of the color to use for the antibiotic.
    bact_name (str): The name of the bacterium.
    melted (bool): A flag that indicates whether the data is in melted format. Default is True.
    """
    
    header=f"""DATASET_STYLE
SEPARATOR COMMA

DATASET_LABEL,{antibiotic}

COLOR,{color}

DATA
"""
    if melted==True:
        criteria=(id_match.drug==antibiotic)&(id_match.resistant=="Resistant")
    else:
        #criteria=id_match[antibiotic]=="1"        
        criteria=id_match[antibiotic]==1
        
    to_label=id_match[criteria]
    #to_label_not=id_match[~criteria]
    
    with open(f"itol_color_{bact_name}_{antibiotic}.txt", "w") as colorfile:
        colorfile.write(header)
        for strain in to_label.panacota_renamed:
            colorfile.write(f"{strain},branch,node,{color},4,normal\n")
            colorfile.write(f"{strain},label,node,#000000,2,bold,{color}\n")    



# Called mutations are loaded
# You need to specify the location of the file with all mutations
# And the file with feature importance analysis result (also a table)

def flatten(inlist):
    #Just list flattening
    outlist=[item for sublist in inlist for item in sublist]
    return outlist

def read_all_mutations(indir,filename):
    #loads and renames a bit, changes type of column
    inpath = os.path.join(indir, filename)
    snp_table=pd.read_csv(inpath, sep="\t")
    snp_table.loc[:,"name/position"]=snp_table["name/position"].astype("str")
    snp_table=snp_table.rename(columns={"name/position":"strain"})
    return(snp_table)

def load_features_to_table(f_dir, snp_table):
    #merges all important features into one table from the all mutations table
    features=[]
    for f_file in os.listdir(f_dir):
        inpath = os.path.join(f_dir, f_file)
        f=pd.read_csv(inpath,sep="\t")
        f_list=f.loc[:,"Mutation name"]
        features.append(f_list) 

    all_f=list(set(flatten(features))) #list of lists become list, which we'll use for column naming
    #we'll also need to save all the antibiotic resistance info
    all_f.extend(["amikacin", "capreomycin", "ethionamide", "kanamycin", "ofloxacin", "streptomycin", "strain"])
    
    f_labels=snp_table.loc[:,all_f]
    return(f_labels)



#snps_dir="/home/sbuyanov/Documents/lab/bacteria_corr/mtb/mtb_called_snps/"
#f_dir="/home/sbuyanov/Documents/lab/bacteria_corr/mtb/feature_importance/without_info_leakage/"

#snp_table=read_all_mutations(snps_dir, "ami_cap_eth_kan_ofl_str_combined_dropped_3.tsv")
#all_f=load_features_to_table(f_dir, snp_table)

#all_f.to_csv("mtb_all_features.tsv", sep="\t")

#snps_dir="/home/sbuyanov/Documents/lab/bacteria_corr/mtb/new_data_feb4/"
#f_dir="/nfs/scistore08/kondrgrp/sbuyanov/muts_in_bact/mtb_strain_data/feature_importance/without_info_leakage/"

#snp_table=pd.read_csv(f"{snps_dir}combined_binary_mutations_non_snp_corrected_0.2_ethionamide.tsv",
#                      sep="\t")