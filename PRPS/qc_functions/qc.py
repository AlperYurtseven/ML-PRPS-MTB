import pandas as pd
from Bio import SeqIO, SearchIO
import os
from tqdm import tqdm
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def make_2000nt_files(strain_file_loc, trunc_fasta_dir):
    """
    Input:
     - strain_file_loc: location of .fasta files
     - trunc_fasta_dir: directory where truncated genome files will be saved

    Output:
     - original file name followed by "_trunc2000.fasta" with truncated genomes.


     Extract first 2000 nucleotides of every contig in every genome file.
     If the sequence is shorter, whole contig is saved

    """

    g_files = os.listdir(strain_file_loc)

    for g_file in tqdm(g_files):
        bact_g = f"{strain_file_loc}/{g_file}"
        with open(f"{trunc_fasta_dir}/{g_file}_trunc2000.fasta", "w") as output_handle:

            for record in SeqIO.parse(bact_g, "fasta"):
                if len(record.seq) <= 2000:
                    SeqIO.write(record, output_handle, "fasta")
                else:
                    record.seq = record.seq[0:2000]
                    SeqIO.write(record, output_handle, "fasta")


def run_blast(fasta_infile, blast_outf, blastdb):
    """
    Uses jupyter magick and runs BLASTN from bash with specified parameters.

    Args:
        fasta_infile (str): Path to FASTA input file.
        blast_outf (str): Output file name for BLASTN results.
        blastdb (str): Path to BLASTN database.

    Returns:
        None.

    Saves the output to a file with the name specified by `blast_outf`.
    """


    !blastn - task
    megablast - db $blastdb - query $fasta_infile - outfmt
    7 - num_alignments
    10 - evalue
    1e-50 - num_threads
    10 > $blast_outf


def run_blast_on_many(trunc_fasta_dir, blast_out_dir, blastdb):
    """
   Runs BLASTN on multiple FASTA input files in a directory.

   Args:
       trunc_fasta_dir (str): Directory containing FASTA input files.
       blast_out_dir (str): Directory where BLASTN output files will be saved.
       blastdb (str): Path to BLASTN database.

   Returns:
       A list of names of all input files that were processed.
   """


os.chdir(trunc_fasta_dir)
non_hits = []
cont_files = os.listdir(trunc_fasta_dir)
for c_file in tqdm(cont_files):
    # infile=strain_file_loc+"/"+c_file
    fasta_infile = trunc_fasta_dir + "/" + c_file
    blast_outf = blast_out_dir + c_file.strip("fasta") + "blast"
    run_blast(fasta_infile, blast_outf, blastdb)
return (cont_files)


def read_blast_out(blast_out_dir):
    """
    Reads BLASTN output files and extracts relevant information.

    Args:
        blast_out_dir (str): Directory containing BLASTN output files.

    Returns:
        A tuple containing:
            - a list of names of empty BLASTN output files
            - a dictionary mapping genome names to the fraction of contigs that did not have a match in the BLASTN database
            - a dictionary mapping genome names to the set of contigs that matched the BLASTN database.
    """
    os.chdir(blast_out_dir)

    wrong_b_cont_d = {}
    empty_blast = []
    mtb_contigs = {}

    outfiles = os.listdir(blast_out_dir)
    for blast_out in tqdm(outfiles):
        no_hits =!grep
        "# 0 hits found" $blast_out | wc - l
        no_hits = int(no_hits[0])

        contigs =!grep
        "hits found" $blast_out | wc - l
        contigs = int(contigs[0])
        if contigs == 0:
            print(f"Empty contig file {blast_out}")
            empty_blast.append(blast_out)
        else:
            wrong_b_cont = no_hits / contigs

            strain = blast_out.split(".fna_prepare")[0]
            wrong_b_cont_d[strain] = wrong_b_cont

            # these are only those which map to mycobacteria
            mtb_contigs[strain] = pd.read_csv(blast_out, comment='#', header=None, sep="\t")[0].unique()

    return (empty_blast, wrong_b_cont_d, mtb_contigs)


def test_num_blast_out():
    """
    Checks that all genome files were successfully processed by BLASTN.
    """
    if empty_blast == 0:
        print("All fasta genome files were blasted")
    if len(empty_blast) + len(wrong_b_cont_d) - len(os.listdir(trunc_fasta_dir)) > 0:
        print("There's a mismatch between files output by blast and number of original files")
    else:
        print("Everything's OK")


def blast_out_dict_to_df(wrong_b_cont_d):
    """
    Converts a dictionary of genome names and their fraction of contigs without a match in the BLASTN database
    to a Pandas DataFrame.

    Args:
        wrong_b_cont_d (dict): A dictionary with genome names as keys and fractions as values.

    Returns:
        A Pandas DataFrame with columns 'strain' (genome name) and 'frac' (fraction of contigs without a match).
        The DataFrame is sorted by descending 'frac' values.
    """
    blast_out_df = pd.DataFrame.from_dict(wrong_b_cont_d, orient="index")
    blast_out_df = blast_out_df.reset_index().rename(columns={0: "frac",
                                                              "index": "strain"})
    # blast_out_df.loc[:,"strain"]=blast_out_df.filename.str.split(".fna_prepare", expand=True)[0]
    blast_out_df = blast_out_df.sort_values(by="frac", ascending=False)
    return blast_out_df