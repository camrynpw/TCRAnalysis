"""
This file should take in an occurrence array for all the HLA alleles
and an occurrence array for all the TCR sequences
and compute where there is overlap between a sequence and an allele
add the counts for number of people with both that allele and that SNE seq, then
compute fisher on that overlap out of total to see which are associated.
It also makes sure to account for NA data in the metadata file and ensures
the filenames are in the same order.
Note this code uses multiprocessing to run this script in parallel.
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from multiprocessing import Pool
import sys

#Read in the metadata file with file information
metadata_file_path = "filepath"
meta_df = pd.read_csv(metadata_file_path, delimiter='\t')
# Account for data that might be missed / NA
missing_allele_files = np.array(meta_df.A.isna())
assert missing_allele_files.sum() == 65
acols = 'A B C DPA1 DPB1 DQA1 DQB1 DRB1 DRB3 DRB4 DRB5'.split()
for col in acols:
    assert all(np.array(meta_df[col].isna())==missing_allele_files)

# Load the file with All Allele bool occurrence lists for all files
allele_file_path = "filepath"
allele_df = pd.read_csv(allele_file_path, header=0)
# Get the allele labels from the first row
alleles = allele_df.columns[1:].tolist()
# Get the allele filenames from the first row as a list of strings
allele_filenames = allele_df.iloc[:, 0].astype(str).tolist()
# Modify filenames as needed to match the metadata file
# Create a new Series in allele_df with the modified filenames
modified_filenames = [f"{filename.replace('.tsv', '')}_imgt.tsv" for filename in allele_filenames]

# Transpose the DataFrame and reset the index
allele_df = allele_df.transpose().reset_index()
# Extract the Boolean occurrence columns
allele_occurrence_array = allele_df.iloc[1:, 1:].values.astype(bool)


# Load the file with the sequence bool occurrence lists for all the files
seq_file_path = "filepath"
seq_df = pd.read_csv(seq_file_path, header=0)

assert seq_df.shape[1] == allele_df.shape[1]
# Get the seq file labels from the first row
seq_filenames = seq_df.columns[1:].astype(str).tolist()
# Assert that the modified filenames in allele_df match the filenames in seq_df
assert seq_filenames == modified_filenames, "File names in allele_df are not in the same format as seq_df."

# Extract the sequence labels from the first column
sequences = seq_df.iloc[:, 0].values
# Extract the Boolean occurrence columns
seq_occurrence_array = seq_df.iloc[:, 1:].values.astype(bool)

# Subset to the columns without missing alleles
seq_occurrence_array = seq_occurrence_array[:,~missing_allele_files]
allele_occurrence_array = allele_occurrence_array[:,~missing_allele_files]
assert seq_occurrence_array.shape[1] == allele_occurrence_array.shape[1]

n_people = seq_occurrence_array.shape[1]
# Assert the missing data is accounted for correctly
assert n_people == 1425 - 65


total_iterations = len(alleles) * len(sequences)
fisher_results = []

# Function to calculate Fisher's exact test and overlaps
def calculate_fisher_overlap(args):
    i, j = args
    n_allele = allele_occurrence_array[i].sum()
    n_seq = seq_occurrence_array[j].sum()
    overlap = (allele_occurrence_array[i] & seq_occurrence_array[j]).sum()

    # Perform Fisher's exact test and store the p-value
    _, p_value = fisher_exact([[overlap, n_allele - overlap], [n_seq - overlap, n_people - n_allele - n_seq + overlap]])

    return (alleles[i], sequences[j], overlap, n_allele, n_seq, p_value)

# Generate all combinations of i and j for multiprocessing
combinations = [(i, j) for i in range(len(alleles)) for j in range(len(sequences))]

# Number of CPU cores to use for multiprocessing
num_cores = 4  # You can adjust this based on your system configuration

# Create a Pool of processes
with Pool(processes=num_cores) as pool:
    # Perform parallel computation
    fisher_results = pool.map(calculate_fisher_overlap, combinations)

# Create a DataFrame from the list of results
results_df = pd.DataFrame(fisher_results, columns=['Allele', 'Sequence', 'Overlap', 'Count_Allele', 'Count_Sequence', 'p_value'])

# Sort the DataFrame by p-value in ascending order
sorted_results_df = results_df.sort_values(by='p_value')

# Save the sorted results to a new CSV file
output_file_path = "filepath"
sorted_results_df.to_csv(output_file_path, index=False)
print("Fisher's test results saved to:", output_file_path)
