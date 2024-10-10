## Script to calculate the background p-values for PyNu triggers

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from multiprocessing import Pool

# Define constants

# Chunk number and GRB name as appears in your .html file
chunk='2'
grbid='grb161210524'
# Directory containing CSV files
input_directory = '/home/hannah.griggs/nu/pynu_tests/o2grbs/results'
results_directory = 'pvals'
# Merged CSV file
merged_csv_path = '{}/matched_{}_{}.csv'.format(input_directory,grbid,'allsky') 
# Z-score parameters
sigma = 500 # Generous +/- 500s window
sigma_index = 500  # Example sigma value for Gaussian index weights
num_target_times = 64000  # Number of different target end times to scan ~box size/30s (~ns inspiral overlap) (estimated for 2 week chunks)

output_csv_path = '{}/{}/top_mod_z_scores_across_target_times_{}.csv'.format(input_directory,
                                                                             results_directory,grbid)
frequency_output_path = '{}/{}/frequency_results_{}.csv'.format(input_directory,
                                                                results_directory,grbid)

# Function to calculate weights based on proximity to the target end time
def calculate_weights(end_times, target_end_time, sigma, index_weights):
    proximity_weights = np.exp(-((end_times - target_end_time) ** 2) / (2 * sigma ** 2))
    return proximity_weights * index_weights

# Function to calculate Gaussian index weights
def calculate_gaussian_index_weights(original_indices, sigma_index):
    return np.exp(-((original_indices - 0) ** 2) / (2 * sigma_index ** 2))

# Calculate the Median Absolute Deviation (MAD)
def calculate_mad(data):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    return mad

# Calculate the modified Z-score for a dataset
def calculate_modified_z_scores(data):
    median = np.median(data)
    mad = calculate_mad(data)
    modified_z_scores = 0.6745 * (data - median) / (mad+0.00001)
    return modified_z_scores

# Function to calculate z-scores
def calculate_z_scores(end_times, weighted_relative_change, mean, std):
    return (weighted_relative_change - mean) / (std+0.1)
# Function to process each target end time
def process_target_end_time(target_end_time, merged_df):
    # Define the range of end times to plot
    end_time_min = target_end_time - 10000
    end_time_max = target_end_time + 10000

    # Filter the DataFrame based on the end time range and conditions
    filtered_df = merged_df[(merged_df['End Time Rounded'] >= end_time_min) & (merged_df['End Time Rounded'] <= end_time_max) &
                            (merged_df['Mass 1_file1'] < 3) & (merged_df['Mass 2_file2'] < 3)]

    # Extract the relevant columns for calculating relative change
    filtered_end_times = filtered_df['End Time Rounded'].values
    ranking_statistic_file1 = filtered_df['Ranking Statistic_file1'].values
    ranking_statistic_file2 = filtered_df['Ranking Statistic_file2'].values

    # Calculate the relative change
    relative_change = ranking_statistic_file2 / ranking_statistic_file1

    # Calculate the weights for the default sigma value
    original_indices = filtered_df['Original Index'].values
    index_weights = calculate_gaussian_index_weights(original_indices, sigma_index)
    weights = calculate_weights(filtered_end_times, target_end_time, sigma, index_weights)
    weighted_relative_change = relative_change * weights

    # Filter relative changes below zero
    #weighted_relative_change = np.clip(weighted_relative_change, 1, None)

    # Calculate the mean and standard deviation of the weighted relative changes
    weighted_relative_change_mean = weighted_relative_change.mean()
    weighted_relative_change_std = weighted_relative_change.std()

    # Calculate the modified Z-scores
    modified_z_scores = calculate_modified_z_scores(weighted_relative_change)

    # Combine end times, Z-scores, and weighted relative changes into a DataFrame
    target_results_df = pd.DataFrame({
        'Target End Time': target_end_time,
        'End Time': filtered_end_times,
        'Weighted Relative Change': weighted_relative_change,
        'Z-score': modified_z_scores
    })

    # Calculate frequency of weighted relative change > 1 for the target end time
    frequency_above_1 = (weighted_relative_change > 1).mean()

    # Get the top 5 Z-scores
    top_5_z_scores_df = target_results_df.nlargest(1, 'Z-score')

    return top_5_z_scores_df, frequency_above_1, target_results_df

# Read the merged CSV file
merged_df = pd.read_csv(merged_csv_path)

# Ensure end times are not printed in scientific notation 
pd.options.display.float_format = '{:.6f}'.format

# Get the range of end times from the merged CSV
end_times_range = merged_df['End Time Rounded'].min(), merged_df['End Time Rounded'].max()
end_times = np.linspace(end_times_range[0], end_times_range[1], num_target_times)

# Add original indices to the merged DataFrame
merged_df['Original Index'] = merged_df.index

# Use multiprocessing to handle target end times
with Pool() as pool:
    results = pool.starmap(process_target_end_time, [(t, merged_df) for t in end_times])

# Collect results
results_df = pd.concat([result[0] for result in results], ignore_index=True)
frequency_target_above_1 = pd.DataFrame({
    'Target End Time': end_times,
    'Frequency Above 1': [result[1] for result in results]
})

# Calculate end time frequencies
end_time_z_scores = pd.concat([result[2] for result in results], ignore_index=True)

# Save the results to a CSV file
results_df.to_csv(output_csv_path, index=False)
print(f"Top Z-scores saved to {output_csv_path}")

# Save the frequency results to a CSV file
frequency_target_above_1.to_csv(frequency_output_path, index=False)
print(f"Frequency of weighted relative change > 1 saved to {frequency_output_path}")

# Calculate the frequency of end_times showing up with the highest or second highest z-score
z_score_ranks = end_time_z_scores.groupby('End Time')['Z-score'].rank(ascending=False)
high_z_score_end_times = end_time_z_scores[z_score_ranks <= 3]

# Frequency of each end time being in the top 2 z-scores across all target end times
frequency_high_z_scores = high_z_score_end_times.groupby('End Time').size().reset_index(name='Frequency')

# Save the frequency results of top z-scores to a CSV file
frequency_high_z_scores_output_path = 'frequency_high_z_scores.csv'
frequency_high_z_scores.to_csv(frequency_high_z_scores_output_path, index=False)
print(f"Frequency of end_times with top Z-scores saved to {frequency_high_z_scores_output_path}")