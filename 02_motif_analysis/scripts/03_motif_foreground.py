####################################################################################
# This script is based on the code:                                                #
# https://github.com/PGB-LIV/Rice_Phospho_Manuscript/blob/main/motif_seqs.py       #       
# creates:                                                                         #  
# 1. foreground = GSB protein sites - 15mers                                       #
# 2. clusterProfiler background = all phosphoproteins in GSB                       #
####################################################################################
import pandas as pd
from Bio import SeqIO
import os
output_location = os.getcwd()
parent_loc = os.path.dirname(output_location)
# This is the all mappings GSB file
phospho_sites = os.path.join(parent_loc, "tgondii_GSB_May2025", "ME49_GSB", "ME49_files", "G2S1B_0.05_protein_pos_single_prot_mapping_ME49.csv")
database = os.path.join(parent_loc, "database", "2024-10-24-decoys-contam-merged_database.fasta")
df = pd.read_csv(phospho_sites)
# Filter for the reference strain (ME49)
df = df[df['Protein'].str.startswith('TGME49')]
# S/T only (matching background)
df = df[df['PTM_residue'].isin(["S", "T"])]
df.to_csv('me49filtered.csv')

#list for all proteins and phosphosites
df['PTM_pos'] = df['Protein_pos'].astype(int)
PTM_pos_list = df['PTM_pos'].tolist()
protein_list = df['Protein'].to_list()
category_list = df['PTM_FLR_category'].to_list()

# Database - full protein sequences
seq_dict = SeqIO.to_dict(SeqIO.parse(database, "fasta"))


###################################################################
#All STY phosphosites +/- 7 - 15mer peptides (motif foreground)   #
#                                                                 # 
###################################################################

bsg_seq_list = []
bsg_seq_protein = []
sg_seq_list = []
sg_seq_protein = []
g_seq_list = []
g_seq_protein = []
bsg_ptm_pos = []
sg_ptm_pos = []
g_ptm_pos = []

# Loop through GSB proteins
for ptm_pos, protein, flr_category in zip(PTM_pos_list,protein_list,category_list):
    original_ptm_pos = ptm_pos
    # check if in database and get seq
    if protein in seq_dict:
        record=seq_dict[protein]
    else:
        print(protein, "not in dict") 
        continue
    # protein seq
    seq_temp=str(record.seq)
    
    # 15mer centered around PTM site
    #remove alanine decoy 
    if seq_temp[ptm_pos - 1] == "A":
        continue
    # if ptm pos +9 greater than seq length (ie not enough aa after ptm pos), add dashes to increase seq length
    if ptm_pos + 9 > len(seq_temp):
        while ptm_pos + 8 > len(seq_temp):
            seq_temp += "_"
    #if ptm pos is less than 9, add dashes before ptm pos (ie want 15mer around ptm pos)
    if ptm_pos < 9:
        seq_temp = ("_" * (8 - ptm_pos)) + seq_temp
        ptm_pos += (8 - ptm_pos)

    # if the ptm is not on STY
    if seq_temp[ptm_pos - 1] != "S" and seq_temp[ptm_pos - 1] != "T":
        print(ptm_pos, protein, flr_category)
        print(seq_temp)
        print(seq_temp[ptm_pos - 8:ptm_pos + 7])
        continue
   #15mer of aa before ptm pos and aa after ptm pos    
   #GSB - list protein, list of seq (15mer) 
    if flr_category == "Bronze":
        bsg_seq_list.append(seq_temp[ptm_pos - 8:ptm_pos + 7])
        bsg_seq_protein.append(protein)
        bsg_ptm_pos.append(original_ptm_pos)
    if flr_category == "Silver":
        bsg_seq_list.append(seq_temp[ptm_pos - 8:ptm_pos + 7])
        bsg_seq_protein.append(protein)
        bsg_ptm_pos.append(original_ptm_pos)
        sg_seq_list.append(seq_temp[ptm_pos - 8:ptm_pos + 7])
        sg_seq_protein.append(protein)
        sg_ptm_pos.append(original_ptm_pos)
    if flr_category == "Gold":
        bsg_seq_list.append(seq_temp[ptm_pos - 8:ptm_pos + 7])
        bsg_seq_protein.append(protein)
        bsg_ptm_pos.append(original_ptm_pos)
        sg_seq_list.append(seq_temp[ptm_pos - 8:ptm_pos + 7])
        sg_seq_protein.append(protein)
        sg_ptm_pos.append(original_ptm_pos)
        g_seq_list.append(seq_temp[ptm_pos - 8:ptm_pos + 7])
        g_seq_protein.append(protein)
        g_ptm_pos.append(original_ptm_pos)
"""
Example of how the padding works
 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19  20 - py index
 M M K L M N A A T S T  D  E  M  K  P  P  T  P  E   R
 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 - actual ptm pos
 ptm occurs at pos 10 on S
 10 - 8 =2 so index 2 is the start
 10 + 7 = 17 (splicing does not include index) so ends at index 16  = 15 mer
  

 not enough aa before ptm
 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  - py index
 N A A T S T D E M K P  P  T  P  E   R
 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 
 ptm at pos 5 on S - ie ptm_pos is less than 9
 8-5 =3 - add 3 dashes at start 
 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17  18- py index
 - - - N A A T S T D  E  M  K  P  P  T  P  E   R
 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
 ptm_pos is now 8 so starts at 0 and ends at 14 index

 not enough aa after ptm 
 0 1 2 3 4 5 6 7 8  - py index
 M K N N A A T S T   
 1 2 3 4 5 6 7 8 9  
 ptm_pos is now 8 on S
 8+9=17
 greater than length of seq
 while seq length is less than 16, keep adding - to the end
 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  - py index
 M K N N A A T S T -  -  -  -  -  - -
 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
 8+7= 15 but it is 14 where ends 
"""

#save as list of sequences - also adding accession
names = ["gsb_motif_seqs.txt", "gs_motif_seqs.txt", "gold_motif_seqs.txt"]
all_seq_list = [bsg_seq_list, sg_seq_list, g_seq_list]
all_proteins = [bsg_seq_protein, sg_seq_protein, g_seq_protein]
ptm_pos =[bsg_ptm_pos, sg_ptm_pos, g_ptm_pos]
for seq, proteins, file_name, ptm_positions in zip(all_seq_list, all_proteins, names, ptm_pos):
    foreground = pd.DataFrame(list(zip(seq, proteins, ptm_positions)), columns=['Sequences','Proteins','PTM_pos'])   
    foreground = foreground.drop_duplicates(["Proteins", "PTM_pos"])
    foreground.to_csv("motif_foreground/"+ file_name, index=False, header=False)


####################################
#All phospho-proteins background   #
#                                  #  
####################################
# Loop through GSB proteins, get UniProt accession, remove duplicates
# gold sites only
df_gold_st = df[df['PTM_FLR_category'] == "Gold"] 
df_gold_st.to_csv("cp_background.csv")
background_protein_list = df_gold_st['Protein'].to_list()

phosphoproteins_df = pd.DataFrame(list(background_protein_list), columns=['Proteins'])
phosphoproteins_df = phosphoproteins_df.drop_duplicates()
phosphoproteins_df.to_csv("cp_background/All_phosphoproteins_background.txt",index=False, header=False)
