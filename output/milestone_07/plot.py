import pandas as pd
from matplotlib import pyplot as plt
import os
import glob
import re

# Customizing the plot for distinct appearance
plt.rcParams["figure.figsize"] = [9.00, 5.50]  # Slightly larger figure size for emphasis
plt.rcParams["figure.autolayout"] = True

# Columns for energy-related data and melting point data
energy_columns = ["iteration", "total_energy", "average_temp", "potential_energy"]
melting_columns = ["cluster_size", "melting_point", "number_of_iterations_till_melting", "added_energy", "atoms_number", "latent_energy"]

directory = './'
file_extension = 'csv'

# Function to extract numeric parts from filenames for sorting
def extract_number_from_filename(filename):
    numbers_in_name = re.findall(r'\d+', filename)
    return int(numbers_in_name[0]) if numbers_in_name else float('inf')

# Function to verify if required columns exist in a file
def check_required_columns(filepath, columns):
    first_row = pd.read_csv(filepath, nrows=1)
    return all(column in first_row.columns for column in columns)

# Change directory and search for energy files recursively
os.chdir(directory)
energy_files = glob.glob(f'*/.{file_extension}', recursive=True)
energy_files.sort(key=extract_number_from_filename)

# Lists for storing legends and heat capacity values
legends = []
heat_capacity_values = []
initial_energies = []

# Loop through each energy file and create plots
for energy_file in energy_files:
    if 'melting_points_vs_cluster_size' in energy_file:
        continue
    if not check_required_columns(energy_file, energy_columns):
        print(f"Skipping {energy_file}: required columns missing.")
        continue

    # Load energy data
    data = pd.read_csv(energy_file, usecols=energy_columns)

    # Plot total energy vs average temperature with a thicker line and new color
    plt.plot(data['average_temp'], data['total_energy'],
             linestyle='-', linewidth=2.5, color='navy')

    legends.append(energy_file.split("/")[0])

    # Calculate heat capacity and add to the list
    delta_energy = data['total_energy'].iloc[-1] - data['total_energy'].iloc[0]
    delta_temp = data['average_temp'].iloc[-1] - data['average_temp'].iloc[0]
    heat_capacity_values.append(delta_energy / delta_temp)
    initial_energies.append(data['total_energy'].iloc[0])

# Update labels, title, and grid styling
plt.xlabel("Average Temp (Kelvin)", fontsize=14, color='darkblue')
plt.ylabel("Total Energy (J)", fontsize=14, color='darkblue')
plt.title("Energy vs Avg Temp - Custom Visualization", fontsize=16, color='black')

# Add a light dashed grid
plt.grid(True, linestyle='--', linewidth=0.7, color='gray')

# Update legend position and font size
plt.legend(legends, loc='upper left', fontsize=9)

# Save the figure and display it
plt.savefig("customized_energy_vs_temp_plot.png", dpi=300)
plt.show()
plt.clf()  # Clear the figure for next plot

# Load melting point data and create the melting point vs cluster size plot
melting_data = pd.read_csv("melting_points_vs_cluster_size.csv", usecols=melting_columns)

# Plot melting point vs cluster size with different marker styles
plt.plot(melting_data['cluster_size'], melting_data['melting_point'],
         linestyle='-', marker='x', color='olive', markersize=8, linewidth=2)

# Customize axes and title
plt.xlabel("Cluster Size", fontsize=13, color='darkgreen')
plt.ylabel("Melting Point (K)", fontsize=13, color='darkgreen')
plt.title("Melting Point vs Cluster Size - Enhanced View", fontsize=15, color='green')

# Add grid with light blue color for distinction
plt.grid(True, linestyle='--', linewidth=0.6, color='lightblue')

# Save the plot with a distinct file name
plt.savefig("melting_point_vs_cluster_size_enhanced.png", dpi=320)
plt.show()
plt.clf()

# Plot heat capacity vs cluster size with updated marker style and line format
plt.plot(melting_data['cluster_size'], heat_capacity_values,
         linestyle='-', marker='d', color='crimson', markersize=9, linewidth=2.2)

# Customize axes, title, and grid
plt.xlabel("Cluster Size", fontsize=13, color='darkred')
plt.ylabel("Heat Capacity (J/K)", fontsize=13, color='darkred')
plt.title("Heat Capacity vs Cluster Size - Detailed", fontsize=15, color='brown')

plt.grid(True, linestyle=':', linewidth=0.6, color='lightgray')

# Save the plot and display it
plt.savefig("heat_capacity_vs_cluster_size_detailed.png", dpi=320)
plt.show()
plt.clf()

# Plot latent heat vs cluster size with new markers and colors
plt.plot(melting_data['cluster_size'], melting_data['latent_energy'],
         linestyle='-', marker='s', color='purple', markersize=8, linewidth=2)

# Update axes and title for the latent heat plot
plt.xlabel("Cluster Size", fontsize=13, color='indigo')
plt.ylabel("Latent Heat (J)", fontsize=13, color='indigo')
plt.title("Latent Heat vs Cluster Size - Refined", fontsize=15, color='darkviolet')

# Change the grid styling for uniqueness
plt.grid(True, linestyle='-', linewidth=0.5, color='lavender')

# Save the final plot and display it
plt.savefig("latent_heat_vs_cluster_size_refined.png", dpi=320)
plt.show()
plt.clf()