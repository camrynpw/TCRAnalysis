"""This file creates occurrence lists for all the different alleles in
the metadata file and saves them to a csv"""

import pandas as pd

# Read the metadata file
metadata_file_path = 'filepath/metadata.csv'
df = pd.read_csv(metadata_file_path, delimiter='\t')

# List of designated columns
designated_columns = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5']

# Initialize an empty DataFrame to store occurrence lists
occurrence_lists_df = pd.DataFrame(index=df.index)

# Placeholder for NA values
na_placeholder = "NA"

# Replace NA values with the placeholder
df[designated_columns] = df[designated_columns].fillna(na_placeholder).astype(str)

# Iterate over each designated column
for col in designated_columns:
    # Split the values on semicolon and create a set of unique values
    unique_values = set(df[col].str.cat(sep=';').split(';'))

    # Remove the placeholder from the unique values set
    unique_values.discard(na_placeholder)

    # Create a dictionary to store the occurrence lists for each unique value
    occurrence_dict = {}
    for value in unique_values:
        occurrence_list = df[col].apply(lambda x: value in x.split(';'))
        occurrence_dict[f"{col}_{value}"] = occurrence_list

    # Concatenate all the occurrence lists for the current designated column
    col_occurrence_df = pd.DataFrame(occurrence_dict)
    occurrence_lists_df = pd.concat([occurrence_lists_df, col_occurrence_df], axis=1)

# Add the "filename" column as the first column by inserting it at the 0th position
occurrence_lists_df = pd.concat([df["filename"], occurrence_lists_df], axis=1)

# Assert to check if the rows under the new "filename" column are in the exact order as filenames from the original metadata file
assert all(occurrence_lists_df["filename"] == df["filename"]), "Rows under the 'filename' column are not in the same order as filenames from the original metadata file."

# Save occurrence lists to a CSV file
output_file_path = 'filepath/All_Allele_occurrence_lists.csv'
occurrence_lists_df.to_csv(output_file_path, index=False)
print("Occurrence lists saved to:", output_file_path)
