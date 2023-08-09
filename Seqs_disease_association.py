"""This script takes in a file that has the occurrence lists for all the sequences
and subsets them by Control and T1D occurrence. Then puts these results into a file that
has the group counts for each sequence. Then runs Fisher on these counts and outputs a csv 
with the pval results"""

import pandas as pd
from scipy.stats import fisher_exact

# Load the CSV files
# The first file needs to contain the bool occurrences for all the sequences for all the files
csv_file_path_1 = 'filepath'
# The second file needs to contain the bool occurrences for the positive group assignment (T1D+) for all the files
csv_file_path_2 = 'filepath'
# The third file needs to contain the bool occurrences for the negative group assignment (T1D-) for all the files
csv_file_path_3 = 'filepath'

df = pd.read_csv(csv_file_path_1, header=0)
df_2 = pd.read_csv(csv_file_path_2)
df_3 = pd.read_csv(csv_file_path_3)

# Extract the sequence labels from the first column
sequences = df.iloc[:, 0].values
# Extract the Boolean occurrence columns
seq_occurrence_array = df.iloc[:, 1:].values.astype(bool)

# Extract the occurrence lists from the second and third files
t1d_occurrence_list = df_2.iloc[:, 0].values.astype(bool)
control_occurrence_list = df_3.iloc[:, 0].values.astype(bool)

# Subset occurrence array based on occurrence lists
t1d_subset_occurrence_array = seq_occurrence_array[:, t1d_occurrence_list]
control_subset_occurrence_array = seq_occurrence_array[:, control_occurrence_list]

# Create column names for subset DataFrames
t1d_column_names = [f"column_{i+1}" for i in range(t1d_subset_occurrence_array.shape[1])]
control_column_names = [f"column_{i+1}" for i in range(control_subset_occurrence_array.shape[1])]

# Create DataFrames for the subsetted occurrence arrays
t1d_subset_df = pd.DataFrame(data=t1d_subset_occurrence_array, columns=t1d_column_names)
control_subset_df = pd.DataFrame(data=control_subset_occurrence_array, columns=control_column_names)

# Insert the Sequence column at the beginning and add count columns
t1d_subset_df.insert(0, "Sequence", sequences)
t1d_subset_df["T1D_Count"] = t1d_subset_df.iloc[:, 1:].sum(axis=1)

control_subset_df.insert(0, "Sequence", sequences)
control_subset_df["NonT1D_Count"] = control_subset_df.iloc[:, 1:].sum(axis=1)

# Merge the DataFrames based on the "Sequence" column
merged_df = pd.merge(t1d_subset_df[["Sequence", "T1D_Count"]], control_subset_df[["Sequence", "NonT1D_Count"]], on="Sequence")

# Run Fisher test
total_T1D = t1d_occurrence_list.sum()
total_control = control_occurrence_list.sum()
assert total_T1D + total_control == 1425

# Create empty lists to store p-values and odds ratios
pvals = []
odds_ratios = []

# Iterate through each row of the DataFrame
for _, row in merged_df.iterrows():
    T1D_count = row["T1D_Count"]
    control_count = row["NonT1D_Count"]

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact([[T1D_count, total_T1D - T1D_count], [control_count, total_control - control_count]], alternative="greater")

    pvals.append(pval)
    odds_ratios.append(oddsratio)

# Add the p-value and odds ratio columns to the DataFrame
merged_df["pval"] = pvals
merged_df["odds_ratio"] = odds_ratios

# Sort the DataFrame by the "pval" column in ascending order
df_sorted = merged_df.sort_values("pval")

# Save the sorted results to a new CSV file
output_file_path = 'filepath'
df_sorted.to_csv(output_file_path, index=False)

print("Fisher's exact test results saved to:", output_file_path)
