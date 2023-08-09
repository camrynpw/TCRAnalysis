"""This file takes in occurrence lists for a variable (such as alleles)
and then occurrence lists for positive and negative group assignments
such as T1D+ and T1D- and subsets the Alleles based on group in order
to calculate association with disease based on Fisher's exact test. 
Results are saved in a csv file with pvals.
"""

import pandas as pd
from scipy.stats import fisher_exact

# First part to subset data frame
# Load the first CSV file with headers with occurrence lists for alleles for all files
csv_file_path_1 = "filepath"
df = pd.read_csv(csv_file_path_1, header=0)

# Get the allele labels from the first row
allele_labels = df.columns[:].tolist()
# Transpose the DataFrame and reset the index
df_1 = df.transpose().reset_index()
# Extract the Boolean occurrence arrray
occurrence_array = df_1.iloc[:, 1:].values.astype(bool)

# Load the t1d CSV file with bool occurrence lists for group assignment (T1D+) for all files
t1d_file_path_2 = "filepath"
df_t1d = pd.read_csv(t1d_file_path_2)
# Extract the occurrence list from the second file
t1d_occurrence_list = df_t1d.iloc[:, 0].values.astype(bool)

# Subset occurrence array based on occurrence list
t1d_subset_occurrence_array = occurrence_array[:, t1d_occurrence_list]

# Create column names for subset_df
t1d_column_names = [f"column_{i+1}" for i in range(t1d_subset_occurrence_array.shape[1])]

# Create a DataFrame for the subsetted occurrence array
t1d_subset_df = pd.DataFrame(data=t1d_subset_occurrence_array, columns=t1d_column_names)
# Add the "Allele" column to the beginning of subset_df
t1d_subset_df.insert(0, "Allele", allele_labels)

# Add a new column "Count" that sums the True values for each row
t1d_subset_df.insert(1, "T1D_Count", t1d_subset_df.iloc[:, 1:].sum(axis=1))

## Do the same for control file
# Load the control CSV file with bool occurrence lists for control group assignment (T1D-) for all files
control_file_path_2 = "filepath"
df_control = pd.read_csv(control_file_path_2)
# Extract the occurrence list from the second file
control_occurrence_list = df_control.iloc[:, 0].values.astype(bool)

# Subset occurrence array based on occurrence list
control_subset_occurrence_array = occurrence_array[:, control_occurrence_list]

# Create column names for subset_df
control_column_names = [f"column_{i+1}" for i in range(control_subset_occurrence_array.shape[1])]

# Create a DataFrame for the subsetted occurrence array
control_subset_df = pd.DataFrame(data=control_subset_occurrence_array, columns=control_column_names)
# Add the "Allele" labels column to the beginning of subset_df
control_subset_df.insert(0, "Allele", allele_labels)

# Add a new column "Count" that sums the True values for each row
control_subset_df.insert(1, "Control_Count", control_subset_df.iloc[:, 1:].sum(axis=1))

# Make sure the dfs have just these two columns
t1d_subset_df = t1d_subset_df[["Allele", "T1D_Count"]]
control_subset_df = control_subset_df[["Allele", "Control_Count"]]

# Merge the DataFrames based on the "Allele" column
merged_df = pd.merge(t1d_subset_df, control_subset_df, on="Allele")


## Run Fisher test
# Define the correct number of subjects
total_T1D = t1d_occurrence_list.sum()
total_control = control_occurrence_list.sum()
assert total_T1D + total_control == 1425

# Create empty lists to store p-values and odds ratios
pvals = []
odds_ratios = []

# Iterate through each row of the DataFrame
for _, row in merged_df.iterrows():
    # Get the counts for T1D and control
    T1D_count = row["T1D_Count"]
    control_count = row["Control_Count"]

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact([[T1D_count, total_T1D - T1D_count], [control_count, total_control - control_count]], alternative="greater")

    # Append the results to the lists
    pvals.append(pval)
    odds_ratios.append(oddsratio)

# Add the p-value and odds ratio columns to the DataFrame
merged_df["pval"] = pvals
merged_df["odds_ratio"] = odds_ratios

# Sort the DataFrame by the "pval" column in ascending order
df_sorted = merged_df.sort_values("pval")

# Save the sorted results to a new CSV file
output_file_path = "filepath"
df_sorted.to_csv(output_file_path, index=False)

print("Fisher's exact test Allele association results saved to:", output_file_path)

