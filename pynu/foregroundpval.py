## Script to calculate the foreground p-values of PyNu triggers

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Z-score parameters
sigma = 500 # Generous +/- 500s window
sigma_index = 500  # Example sigma value for Gaussian index weights
num_target_times = 64000  # Number of different target end times to scan ~box size/30s (~ns inspiral overlap) 

# GRB information
grbid='grb161210524'
center_end_time = 1165408451 # GRB T0
trigger_timewindow = 10  # Generous physically-motivated time delay betwen GW and GRB

# Load your data
# Directory containing CSV files
input_directory = '/home/hannah.griggs/nu/pynu_tests/o2grbs/results/pvals'
results_df = pd.read_csv('{}/top_mod_z_scores_across_target_times_{}.csv'.format(input_directory,grbid))
timewindow=sigma*2

# Function to calculate the average Z-score for each 'End Time'
def calculate_average_z_scores(results_df):
    average_z_scores = results_df.groupby('End Time')['Z-score'].mean()
    return average_z_scores

# Function to calculate the average Z-score for each 'End Time' where weighted_relative_change > 1
def calculate_filtered_average_z_scores(results_df):
    filtered_df = results_df[results_df['Weighted Relative Change'] > 1]
    average_z_scores = filtered_df.groupby('End Time')['Z-score'].mean()
    return average_z_scores

# Function to calculate the frequency of each 'End Time' being a top Z-score
def calculate_top_z_score_frequencies(results_df):
    filtered_df = results_df[results_df['Weighted Relative Change'] > 1.1]
    top_z_scores = filtered_df.groupby('Target End Time')['Z-score'].idxmax()
    top_z_scores_df = filtered_df.loc[top_z_scores]
    top_z_scores_frequencies = top_z_scores_df['End Time'].value_counts()
    return top_z_scores_frequencies

# Function to adjust average Z-scores by multiplying with frequency squared
def adjust_z_scores_by_frequency_squared(average_z_scores, top_z_scores_frequencies):
    adjusted_z_scores = average_z_scores.copy()
    for end_time, frequency in top_z_scores_frequencies.items():
        if end_time in adjusted_z_scores.index:
            adjusted_z_scores[end_time] *= 1#frequency ** 1
    return adjusted_z_scores

# Function to calculate the chance that another end time has a top Z-score >= signal's top Z-score
def calculate_chance_of_higher_or_equal_z_score(results_df, top_signal_z_score, center_end_time):
    #print(len(results_df['Target End Time']))
    # Exclude the signal end time
    non_signal_df = results_df[results_df['Target End Time'] != center_end_time]
    #print(len(non_signal_df['Target End Time']))
    # Count how many non-signal end times have a top Z-score >= the signal's top Z-score
    higher_or_equal_z_scores = non_signal_df[non_signal_df['Z-score'] >= top_signal_z_score]
    num_higher_or_equal = len(higher_or_equal_z_scores)

    # Calculate the total number of unique non-signal end times
    total_non_signal_end_times = non_signal_df['Target End Time'].nunique()

    # Calculate the probability
    probability = num_higher_or_equal / total_non_signal_end_times if total_non_signal_end_times > 0 else 0

    return probability, num_higher_or_equal, total_non_signal_end_times

# Function to find the top Z-score for the signal end time
def find_top_signal_z_score(results_df, center_end_time, time_window):
    # Find rows where 'End Time' is within the time window around center_end_time
    signal_df = results_df[(results_df['End Time'] >= center_end_time - time_window) &
                           (results_df['End Time'] <= center_end_time + time_window)]

    # If no rows are found, return NaN for both Z-score and End Time
    if signal_df.empty:
        print(f"No trigger found within +/- {time_window} seconds of the center_end_time: {center_end_time}")
        return np.nan, np.nan

    # Find the row with the maximum Z-score
    top_signal_row = signal_df.loc[signal_df['Z-score'].idxmax()]

    # Extract the top Z-score and the corresponding End Time
    top_signal_z_score = top_signal_row['Z-score']
    new_center_end_time = top_signal_row['End Time']

    return top_signal_z_score, new_center_end_time

# Function to plot adjusted Z-scores and frequencies on separate plots
def plot_adjusted_z_scores_and_frequencies(adjusted_z_scores, top_z_scores_frequencies, center_end_time):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

    # Plot adjusted Z-scores
    centered_x = adjusted_z_scores.index - center_end_time
    ax1.plot(centered_x, adjusted_z_scores.values, 'bo', label='Adjusted Z-Scores')
    ax1.scatter([0], adjusted_z_scores.loc[center_end_time], color='blue', s=100, edgecolor='black', zorder=5, label='Highlighted Physical Signal')
    ax1.axvline(0, color='black', linestyle='--', lw=1, label=f'Center: {center_end_time}')
    ax1.set_ylabel('Adjusted Z-score')
    ax1.set_title('Adjusted Z-scores for Each End Time')
    ax1.legend()
    ax1.grid(True)

    # Plot frequencies
    centered_frequencies_x = top_z_scores_frequencies.index - center_end_time
    if not centered_frequencies_x.empty:
        ax2.scatter(centered_frequencies_x, top_z_scores_frequencies.values,color='blue', alpha=0.7, label='Frequency of Top Z-Score')
        ax2.set_xlabel(f'End Time Relative to {center_end_time}')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Frequency of Each End Time as Top Z-Score (Weighted Relative Change > 1)')
        ax2.axvline(0, color='black', linestyle='--', lw=1, label=f'Center: {center_end_time}')
        ax2.legend()
        ax2.grid(True)
    else:
        ax2.text(0.5, 0.5, 'No data available for frequencies', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, fontsize=12, color='red')

    plt.tight_layout()
    plt.show()

# Calculate the average Z-score for each 'End Time' where weighted_relative_change > 1
filtered_average_z_scores = calculate_filtered_average_z_scores(results_df)

# Calculate the frequencies of each 'End Time' being the top Z-score where weighted_relative_change > 1
top_z_scores_frequencies = calculate_top_z_score_frequencies(results_df)

# Adjust the average Z-scores by multiplying with the frequency squared
adjusted_z_scores = adjust_z_scores_by_frequency_squared(filtered_average_z_scores, top_z_scores_frequencies)

# Plot the adjusted Z-scores and frequencies
#plot_adjusted_z_scores_and_frequencies(adjusted_z_scores, top_z_scores_frequencies, center_end_time)

# Find the top Z-score for the signal end time
top_signal_z_score, new_center_end_time = find_top_signal_z_score(results_df, center_end_time, trigger_timewindow)

# Calculate the chance of a higher or equal Z-score for non-signal end times
probability, num_higher_or_equal, total_non_signal_end_times = calculate_chance_of_higher_or_equal_z_score(results_df, top_signal_z_score, new_center_end_time)

print(f"Top Z-score for signal end time {center_end_time}: {top_signal_z_score} at {new_center_end_time}")
print(f"Number of non-signal end times with Z-score >= {top_signal_z_score}: {num_higher_or_equal}")
#print(f"Total number of unique non-signal end times: {total_non_signal_end_times}")
print(f"Probability of another end time having a Z-score >= {top_signal_z_score}: {probability:.4f}")

