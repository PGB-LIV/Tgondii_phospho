######################################################
# This script is for generating 15mer background     # 
#                                                    #   
######################################################
import sys
import pandas as pd
from Bio import SeqIO
import os
import re

# Database used in analysis
database = "2024-10-24-decoys-contam-merged_database.fasta"
# into dictionary
seq_dict = SeqIO.to_dict(SeqIO.parse(database, "fasta"))

# input and outputs for each exp
fdr_input = sys.argv[1]
background_output= sys.argv[2]

overall_background = []

# For each exp...
# FDR_output.csv from mzidFLR
all_sites = fdr_input + "/FDR_0.01/FDR_output.csv"
df = pd.read_csv(all_sites)
# <1% FDR
df = df.loc[df['FDR'] <= float(0.01)]
# Drop duplicates on protein - peptide
df = df.drop_duplicates(subset=["Protein","Peptide"], keep="first").reset_index(drop=True)
protein_list = df['Protein'].to_list()
peptide_list = df['Peptide'].to_list()
# Loop through protein/peptide
for peptide, protein in zip(peptide_list, protein_list):
    # if peptide maps to multiple proteins, keep first only
    all_proteins = protein.split(":")
    # index of S and T in peptide
    pep_index = [st.start() for st in re.finditer("[ST]", peptide)]
    
    # ME49 only
    me49_proteins = [prot for prot in all_proteins 
                  if "TGME49" in prot and "rev_" not in prot]
    # single mapping: take first ME49 protein only, skip if none
    if len(me49_proteins) == 0:
        continue
    protein = me49_proteins[0]

    # Look up full protein seq in dict
    record = seq_dict[protein]
    seq_temp = str(record.seq)
  
    # If there is at least 1 S|T
    if len(pep_index) > 0:
        # find index of peptide in protein
        prot_index = [peptide_pos.start() for peptide_pos in re.finditer(peptide, seq_temp)]
        if len(prot_index) == 0:
            print("mapping error for: ", peptide, " in", protein)
        # for each instance of the peptide in the protein (I guess this wont occur very often)
        for index in prot_index:
            # create new list with index of S|T in the protein (index of S|T in peptide + index of peptide in protein)
            aa_index_in_prot = [st_index + index for st_index in pep_index]
            # For each index of S|T in the protein (from identified peptides)
            for aa in aa_index_in_prot:
                true_aa_pos = aa
                #if len(prot_index)>1:
                    #print(aa,index,prot_index, peptide, seq_temp) 
                # clean protein seq
                seq_temp = str(record.seq)
                if seq_temp[aa] != "S" and seq_temp[aa] != "T":
                    print("15mer fail from", peptide, "/", seq_temp[aa])
                    continue
                # needs to be at least 7 aa after S|T, so add "_" to end of seq
                if aa + 8 > len(seq_temp):
                    while aa + 8 > len(seq_temp):
                        seq_temp += "_"
                # needs to be at least 7 aa before S|T, so add "_" before seq
                if aa < 7:
                    seq_temp = ("_" * (7 - aa)) + seq_temp
                    aa += (7 - aa) # update index of S|T in protein since the dashes will shift the position
                    if seq_temp[aa] != "S" and seq_temp[aa] != "T":
                        print("15mer fail 2", seq_temp[aa - 7:aa + 8])
                        continue
                # string of 15mer and protein
                protein_str = record.id + ":" + seq_temp[aa - 7:aa  + 8] + ":" + str(true_aa_pos + 1) #py is 0 based   
                if protein_str not in overall_background:
                    overall_background.append(protein_str)    
    
print(len(overall_background))
df2=pd.DataFrame(list(zip(overall_background)),columns=["Proteins"])
# Split protein and 15mer into columns   
df2[["Proteins","Sequence", "PTM_protein_position"]]=df2["Proteins"].str.split(":", n=2,expand=True)                               

df2 = df2[["Sequence","Proteins", "PTM_protein_position"]]
df2.to_csv(background_output, index=False)
    