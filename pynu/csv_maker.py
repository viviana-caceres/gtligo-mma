import os
import pandas as pd
from bs4 import BeautifulSoup
os.environ['OPENBLAS_NUM_THREADS'] = '1'

## Initialize variables

# Directory containing HTML files
input_directory = '/home/hannah.griggs/nu/pynu_tests/o2grbs/results/'

# Chunk number and GRB name as appears in your .html file
chunk='2'
grbid='grb161210524'

def extract_data_from_html(file_path):
    with open(file_path, 'r') as f:
        html_content = f.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    # Find the <script> tag that contains the data
    script_tags = soup.find_all('script')
    data_string = ""

    for script in script_tags:
        if script.string and 'data.addRows' in script.string:
            data_string = script.string
            break

    ## Extract the data string
    if data_string:
        start_index = data_string.index('addRows(') + len('addRows(')
        end_index = data_string.index(');', start_index)
        data_rows = data_string[start_index:end_index].strip()
    
        # Convert the string representation of rows into a list of lists
        rows = eval(data_rows)
    
        # Extract relevant entries
        exc_ifar = [row[0] for row in rows]  # Exc. IFAR (YR)
        inc_ifar = [row[1] for row in rows]  # Inc. IFAR (YR)
        exc_fap = [row[2] for row in rows]    # Exc. FAP
        inc_fap = [row[3] for row in rows]    # Inc. FAP
        ranking_statistics = [row[4] for row in rows]  # Ranking Statistic
        end_times = [row[5] for row in rows]           # End Time
        m1 = [row[7] for row in rows]   # Mass 1
        m2 = [row[8] for row in rows]   # Mass 2
        h1chi2 = [row[12] for row in rows]   # H1 Red Chi2
        l1chi2 = [row[14] for row in rows]   # L1 Red Chi2
    
        # Create a DataFrame
        df = pd.DataFrame({
            'End Time': end_times,
            'Exc. IFAR (YR)': exc_ifar,
            'Inc. IFAR (YR)': inc_ifar,
            'Exc. FAP': exc_fap,
            'Inc. FAP': inc_fap,
            'Ranking Statistic': ranking_statistics,
            'Mass 1': m1,
            'Mass 2': m2,
            'H1 Red. Chisq': h1chi2,
            'L1 Red. Chisq': l1chi2
        })

    return df

## Create a CSV file containing the relevant results from PyCBC's output
# Process each HTML file in the directory
for filename in os.listdir(input_directory):
    if filename.endswith('.html'):
        csv_filename = filename.replace('.html', '.csv')
        csv_path = os.path.join(input_directory, csv_filename)

        # Check if the CSV file already exists
        if not os.path.exists(csv_path):
            file_path = os.path.join(input_directory, filename)
            df = extract_data_from_html(file_path)
            df.to_csv(csv_path, index=False)

print("Processing complete. CSV files created or already existing files were skipped.")

## Create a merged CSV file with the results from the all-sky search.
csv_file1 = input_directory+'outputallskychunk'+chunk+'_FG.csv'
csv_file2 = input_directory+'output'+grbid+'_FG.csv'

# Read the CSV files into DataFrames
df1 = pd.read_csv(csv_file1)
df2 = pd.read_csv(csv_file2)
# Round 'End Time' values to three decimal places
df1['End Time Rounded'] = df1['End Time'].round(1)
df2['End Time Rounded'] = df2['End Time'].round(1)

# Find common 'End Time Rounded' values
common_end_times = set(df1['End Time Rounded']).intersection(set(df2['End Time Rounded']))

# Filter the DataFrames to keep only the rows with common 'End Time Rounded' values
df1_common = df1[df1['End Time Rounded'].isin(common_end_times)]
df2_common = df2[df2['End Time Rounded'].isin(common_end_times)]

# Merge the DataFrames on 'End Time Rounded'
merged_df = pd.merge(df1_common, df2_common, on='End Time Rounded', suffixes=('_file1', '_file2'))

# Save the merged DataFrame to a new CSV file
output_csv_path = '{}/matched_{}_{}.csv'.format(input_directory,grbid,'allsky')
merged_df.to_csv(output_csv_path, index=False)

print("Comparison complete. Merged CSV file created.")