import pandas as pd
from matplotlib import pyplot as plt
import os

# Define the path to the CSV file
csv_file_path = 'energies_vs_time_step.csv'

# Check if the CSV file exists
if not os.path.exists(csv_file_path):
    raise FileNotFoundError(f"The file {csv_file_path} does not exist.")

# Configure plot settings
plt.rcParams["figure.figsize"] = [8.00, 4.00]  # Set figure size
plt.rcParams["figure.autolayout"] = True  # Ensure layout fits into the figure area

# Specify the columns to read from the CSV file
columns = ["time_step", "total_energy", "kinetic_energy", "potential_energy"]

# Load the data from the CSV file
df = pd.read_csv(csv_file_path, usecols=columns)

# Create the plot
plt.figure()

# Plot total energy with updated style and color
plt.plot(df.time_step, df.total_energy,
         label='Total Energy',
         color='teal',
         linestyle='--',
         marker='o',
         markersize=6,
         linewidth=2)

# Add labels and title
plt.xlabel("Time Step", fontsize=12)
plt.ylabel("Total Energy", fontsize=12)
plt.title("Total Energy vs Time Step", fontsize=14)

# Add grid for better readability
plt.grid(True, linestyle='--', alpha=0.7)

# Add legend to the plot
plt.legend(loc='best')

# Save the plot as a PNG file
plt.savefig("total_energy_vs_time_step_updated.png", dpi=300)

# Display the plot
plt.show()
