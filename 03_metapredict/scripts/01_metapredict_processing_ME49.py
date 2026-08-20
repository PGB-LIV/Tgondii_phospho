########################
# Disorder analysis    #
# Tgondii PHOSPHO      #  
########################
import pandas as pd
import Bio.SeqIO as SeqIO
import os
############################################################
# 1. Get max seq lenghth                                   #
# Needed for reading in df correctly                       # 
# with column header for index of score                    #
# The csv has no header and the rows are different lengths # 
# col 1 = protein                                          #    
# col 2= seq                                               # 
# col 3-n = disorder scores                                #
############################################################

# create results dir
os.makedirs("metapredict_disorder_processed_ME49", exist_ok=True)
# get max seq length
# read in first 2 columns
seq_len = pd.read_csv("ME49tgondii_metapredict_v3_240726.csv", usecols=[0,1], header=None)
seq_len.columns =["protein", "seq"]
# rm spaces
seq_len['seq'] = seq_len['seq'].str.strip()
# col with seq length
seq_len['seq_len'] = seq_len['seq'].str.len()
seq_len.to_csv("metapredict_disorder_processed_ME49/seqlen.csv")
# Longest seq (ie column header needs to go up to this number)
max_seq_length = seq_len['seq_len'].max()
print("Max sequence length:", max_seq_length)

# Column names
col_list = ["Protein", "Sequence"]
# col names of index pos
pos_list = list(range(1, max_seq_length + 1))
col_list = col_list + pos_list

###############################
# 2. Now read in full df with #
# column names as above       #   
###############################

# Read in data with specificied col names
df = pd.read_csv("ME49tgondii_metapredict_v3_240726.csv", header=None, names=col_list) # From metapredict (on cluster), need to gives names due to different lengths of rows
# remove spaces from seq column
df['Sequence'] = df['Sequence'].str.strip()
# separate protein id
df[["Protein","Description"]] = df["Protein"].str.split(" ", n=1,expand=True)
df = df.drop('Description', axis=1)
# this is ME49 only from ME49 only DB (no contam or decoy) --> no need to remove these in this script
# Wide -> long
df = pd.melt(df, id_vars=["Protein","Sequence"])
df.to_csv("metapredict_disorder_processed_ME49/melt.csv")
# remove blanks (ie proteins with shorter seq)
df = df.dropna(subset=["value"])  
df = df.rename(columns={"variable":"Protein_pos", "value":"disorder"})

# get aa corresponding to row
df["PTM_residue"] = df.apply(
    lambda row: row["Sequence"][row["Protein_pos"] - 1] if pd.notnull(row["Sequence"]) else None,
    axis=1
)
# ST only
df = df[(df["PTM_residue"]=="S")|(df["PTM_residue"]=="T")]

############################
# 3. merge with GSB file   #
#                          # 
############################ 

output_location = os.getcwd()
parent_loc = os.path.dirname(output_location)         
# add col to say if aa is a phosphosite in gsb
gsb = pd.read_csv(os.path.join(parent_loc, "tgondii_GSB_May2025", "ME49_GSB", "ME49_files", "G2S1B_0.05_protein_pos_single_prot_mapping_ME49.csv"))
# S or T only
gsb = gsb[(gsb["PTM_residue"]=="S")|(gsb["PTM_residue"]=="T")]
gsb = gsb[["Protein", "Protein_pos", "PTM_residue", "PTM_FLR_category"]]
gsb["Phospho"] = True
df = df.merge(gsb, how="left", on=["Protein", "Protein_pos", "PTM_residue"])
df["Phospho"] = df["Phospho"].fillna(False)
df=df[["Protein", "Protein_pos", "disorder", "PTM_residue", "PTM_FLR_category", "Phospho"]]
df.to_csv("metapredict_disorder_processed_ME49/disorder_processed_tgondii.csv", index=False)

