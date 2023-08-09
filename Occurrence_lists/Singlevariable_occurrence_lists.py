import pandas as pd

# Read the metadata file, (this file is actually tsv)
metadata_file_path = 'filepath/metadata.tsv'
metadata_df = pd.read_csv(metadata_file_path, delimiter='\t')

# Making lists from metadata file. Make sure in the right order

T1D_occurrence_list = [status == 'T1D' for status in metadata_df['diabetes_status']]
Control_occurrence_list = [status == 'CTRL' for status in metadata_df['diabetes_status']]
Family_occurrence_list = [(status == 'FDR') or (status == 'SDR') for status in metadata_df['diabetes_status']]
nonT1D_occurrence_list = [(status == 'FDR') or (status == 'SDR') or (status == 'CTRL') for status in metadata_df['diabetes_status']]

# Create DataFrames for the occurrence lists and save them
T1D_df = pd.DataFrame({'Occurrence': T1D_occurrence_list})
T1D_output_file_path = 'filepath/T1D_occurrence_list.csv'
T1D_df.to_csv(T1D_output_file_path, index=False)

Control_df = pd.DataFrame({'Occurrence': Control_occurrence_list})
Control_output_file_path = 'filepath/Control_occurrence_list.csv'
Control_df.to_csv(Control_output_file_path, index=False)

Family_df = pd.DataFrame({'Occurrence': Family_occurrence_list})
Family_output_file_path = 'filepath/Family_occurrence_list.csv'
Family_df.to_csv(Family_output_file_path, index=False)

nonT1D_df = pd.DataFrame({'Occurrence': nonT1D_occurrence_list})
nonT1D_output_file_path = 'filepath/nonT1D_occurrence_list.csv'
nonT1D_df.to_csv(nonT1D_output_file_path, index=False)

print("Occurrence lists saved as CSV files.")
