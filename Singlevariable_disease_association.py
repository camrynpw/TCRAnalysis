import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

# First part to subset data frame
# Load the first CSV file with headers that has a single occurrence list for a variable (example: 1 allele)
csv_file_path_1 = "filepath"
df = pd.read_csv(csv_file_path_1, header=0)
# Transpose the DataFrame and reset the index
df_1 = df.transpose().reset_index()
# Extract the Boolean occurrence arrray
occurrence_array = df_1.iloc[:, 1:].values.astype(bool)
antigen_occs = np.array(occurrence_array)

# Load the t1d CSV file with bool occurrence lists for control group assignment (T1D+) for all files 
t1d_file_path = "filepath"
df_t1d = pd.read_csv(t1d_file_path)
t1d_occurrence_list = df_t1d.iloc[:, 0].values.astype(bool)
t1d_occs = np.array(t1d_occurrence_list)

# Subset occurrence array based on occurrence list
t1d_subset_occs = antigen_occs & t1d_occs
t1d_count = t1d_subset_occs.sum()

## Do the same for control file
# Load the control CSV file with bool occurrence lists for control group assignment (T1D-) for all files
control_file_path = "filepath"
df_control = pd.read_csv(control_file_path)
control_occurrence_list = df_control.iloc[:, 0].values.astype(bool)
control_occs = np.array(control_occurrence_list)

# Subset occurrence array based on occurrence list
control_subset_occs = antigen_occs & control_occs
control_count = control_subset_occs.sum()


## Run Fisher test
# Define the total number of subjects in T1D and control groups
total_T1D = t1d_occs.sum()
total_control = control_occs.sum()

# Perform Fisher's exact test
_, pval = fisher_exact([[t1d_count, total_T1D - t1d_count], [control_count, total_control - control_count]], alternative="greater")


print("Variable association with positive group assignment (T1D+): ", pval)

