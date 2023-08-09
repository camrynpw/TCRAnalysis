''' 
This script takes in a list of TCR sequences of interest (ones to be analyzed)
and figures out which files/subjects it occurs in to create occurrence lists
of bools for every file.
And it saves a CSV file with all that occurrence information in it. The format of that file
is: one row per TCR; one column per subject (ie 1425); True or False to indicate if TCR occurs
in that subject
Note that the column names are modified to be the filenames and makes sure the files are processed
in the same order so the lists will line up for further analysis.
'''

import pandas as pd
import os
import numpy as np

# Read the unique importance sequences from the file. Must have a column named "Sequence"
sequences_file_path = "filepath/sequences.csv"
sequences_df = pd.read_csv(sequences_file_path)

# Define the folder path containing the cohort files
cohort_folder_path = "filepath/cohort_files"

tcrs = sequences_df['Sequence'].tolist()

num_tcrs = len(tcrs)
num_files = 1425

# create a dictionary that maps from tcr to integer index
tcr2int = {t: i for i, t in enumerate(tcrs)}  # dictionary comprehension

# Creates the occurrence matrix of falses
all_occs = np.zeros((num_tcrs, num_files), dtype=bool)

# Sort the filenames in the exact order found in the metadata file
# Step 1: Read the CSV file into a DataFrame
csv_file_path = "filepath/metadata.csv"
metadata_df = pd.read_csv(csv_file_path, delimiter='\t')

# Step 2: Get the filenames from a specific column in the DataFrame and modify them to match
filename_list = metadata_df['filename'].tolist()
#Fix the filename adjustment as needed
filenames = [f"{os.path.splitext(filename)[0]}_imgt.tsv" for filename in filename_list]


for filenum, filename in enumerate(filenames):
    file_path = os.path.join(cohort_folder_path, filename)
    
    # Read the TCRs in the repertoire from the file
    file_df = pd.read_csv(file_path, delimiter='\t')
    file_tcrs = file_df['combined_seq'].tolist()

    # Go through and record the occurrence for each TCR
    for tcr in file_tcrs:
        if tcr in tcr2int:
            # Record the occurrence
            all_occs[tcr2int[tcr], filenum] = True

# Convert the occurrence matrix to a DataFrame
occurrence_df = pd.DataFrame(all_occs, columns=filenames, index=tcrs)

assert all(occurrence_df.columns == filenames), "Columns in occurrence_df do not match filename_list order."

# Save the occurrence DataFrame as a CSV file
output_file_path = "filepath/allseq_occs.csv"
occurrence_df.to_csv(output_file_path, index=True)

print("Sequence occurrence array saved to:", output_file_path)
