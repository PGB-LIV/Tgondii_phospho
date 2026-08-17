####################################
# 15mer background                 #
# script 2                         #
# combine 15mers for each exp      #  
# and remove duplicates            #
####################################
import pandas as pd
import glob
#https://stackoverflow.com/questions/77238010/combine-multiple-csv-files-using-merge-while-retaining-filename-information

# All csv files in directory
background_exp_files = glob.glob('background_per_file/*.csv')

# concatenate - reads file and stores in list to then combine
df = pd.concat([pd.read_csv(file) for file in background_exp_files], ignore_index=True)

df.to_csv("combined_all.csv",index=False)
print(len(df))
# Remove duplicates on protein and pos
df=df.drop_duplicates(subset=["Proteins", "PTM_protein_position"])
print(len(df))
df.to_csv("combined_final.csv",index=False)

# Final txt file
df.to_csv("background_combined/motif_ST_background.txt",index=False, header=False)

