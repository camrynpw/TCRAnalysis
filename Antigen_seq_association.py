import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

# Read the metadata file, remember this file is actually tsv
metadata_file_path = 'filepath/metadata.csv'
metadata_df = pd.read_csv(metadata_file_path, delimiter='\t')

#Make a list for each antigen column
GAD_list = [(status == 1.0) for status in metadata_df['GAD_Status']]
GAD_occs = np.array(GAD_list)

IA2_list = [(status == 1.0) for status in metadata_df['IA2_Status']]
IA2_occs = np.array(IA2_list)

ZnT8_list = [(status == 1.0) for status in metadata_df['ZnT8_Status']]
ZnT8_occs = np.array(ZnT8_list)

a_ab_list = [(status == 3.0) for status in metadata_df['a_ab']]
a_ab_occs = np.array(a_ab_list)

# Create DataFrames for the occurrence lists and save them
GAD_df = pd.DataFrame({'Occurrence': GAD_list})
GAD_file_path = "filepath/GAD_occurrence_list.csv"
GAD_df.to_csv(GAD_file_path, index=False)

IA2_df = pd.DataFrame({'Occurrence': IA2_list})
IA2_file_path = "filepath/IA2_occurrence_list.csv"
IA2_df.to_csv(IA2_file_path, index=False)

ZnT8_df = pd.DataFrame({'Occurrence': ZnT8_list})
ZnT8_file_path = "filepath/ZnT8_occurrence_list.csv"
ZnT8_df.to_csv(ZnT8_file_path, index=False)

a_ab_df = pd.DataFrame({'Occurrence': a_ab_list})
a_ab_file_path = "filepath/a_ab_occurrence_list.csv"
a_ab_df.to_csv(a_ab_file_path, index=False)


all_three_occs = GAD_occs & IA2_occs & ZnT8_occs
assert all_three_occs.sum() == a_ab_occs.sum()

# Now the code will be comparing antigen occurrences with the TCR sequences
# Read in the sequence bool occurrences for all files
seq_file_path = 'filepath/allseq_occs.csv'
seq_df = pd.read_csv(seq_file_path, header=0)

seq_df = seq_df.transpose()
seq_df.columns = seq_df.iloc[0]
seq_df = seq_df.drop('Unnamed: 0')

n_people = 1425

# Create an empty list to store the results
results = []

# Run through all the sequences in the seq file and find associations for them
for sequence in seq_df.columns:

    seq_occs = seq_df[sequence]
    seq_occs = seq_occs.astype(bool)
    seq_occs = np.array(seq_occs)

    assert a_ab_occs.shape == seq_occs.shape

    TCR_within_antigen_occs = seq_occs & a_ab_occs

    n_antigen = a_ab_occs.sum()
    n_seq = seq_occs.sum()
    overlap = TCR_within_antigen_occs.sum()

    # Perform Fisher's exact test and store the p-value
    _, p_value = fisher_exact([[overlap, n_antigen - overlap], [n_seq - overlap, n_people - n_antigen - n_seq + overlap]])

    # Save the result if p-value < 0.05
    if p_value < 0.05:
        results.append({'Sequence': sequence, 'p_value': p_value, 'Overlap': overlap})

# Create a DataFrame from the list of results
result_df = pd.DataFrame(results)

# Sort the DataFrame by p-value in ascending order
result_df = result_df.sort_values(by='p_value')

# Save the result_df to a CSV file
result_file_path = "filepath/Triple_antigen_seq_association.csv"
result_df.to_csv(result_file_path, index=False)

print("Results saved to: ", result_file_path)
