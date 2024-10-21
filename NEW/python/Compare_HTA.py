import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Step 1: Read CSV files
df1 = pd.read_csv('event_info.csv')
df2 = pd.read_csv('event_info_HTA.csv')
df3 = pd.read_csv('forJunewoo_ITA_SR_events.csv')

# Step 2: Create composite keys from 'experiment', 'run', and 'event'
df1['composite_key'] = df1.apply(lambda row: (row['experiment'], row['run'], row['event']), axis=1)
df2['composite_key'] = df2.apply(lambda row: (row['experiment'], row['run'], row['event']), axis=1)
df3['composite_key'] = df3.apply(lambda row: (row['experiment'], row['run'], row['event']), axis=1)

# Step 3: Get unique composite keys
keys1 = set(df1['composite_key'])
keys2 = set(df2['composite_key'])
keys3 = set(df3['composite_key'])

# Step 4: Calculate the number of unique and shared composite keys
#only_in_file1 = len(keys1 - keys2)  # Unique to file1
#only_in_file2 = len(keys2 - keys1)  # Unique to file2
#shared_keys = len(keys1 & keys2)    # Shared in both files

# Step 5: Draw the Venn diagram
#venn2(subsets=(only_in_file1, only_in_file2, shared_keys),
#      set_labels=(r'$B \rightarrow X_{s} \nu \bar{\nu}$', r'$B^{+} \rightarrow K^{+} \nu \bar{\nu}$ ITA'))


# Step 5: Create a Venn diagram
plt.figure(figsize=(8, 8))
venn3([keys1, keys2, keys3], ('Group 1', 'Group 2', 'Group 3'))
# Display the Venn diagram
plt.show()
