import pandas as pd
from matplotlib import pyplot as plt

atoms_number = 923

# Setting up plot parameters with a new style for unique visualization
plt.rcParams["figure.figsize"] = [9.00, 5.50]  # Adjusted plot size for clarity
plt.rcParams["figure.autolayout"] = True

# First Plot: Demonstrating Energy Conservation Time vs Core Count
data_columns = ["cores", "real", "user", "sys"]
df = pd.read_csv(f"{atoms_number}/equilibrium_on_923_atoms.csv", usecols=data_columns)

# Plotting User Time vs Cores with distinct style and color
plt.plot(df['cores'], df['user'], linestyle='-', color='navy', linewidth=2.2)
plt.xlabel("Number of Cores", fontsize=13, fontweight='medium', color='navy')
plt.ylabel("Energy Conservation Time (fs)", fontsize=13, fontweight='medium', color='navy')
plt.title("Energy Conservation Time vs Core Count", fontsize=15, fontweight='bold', color='darkred')

# Customizing the grid for uniqueness
plt.grid(True, linestyle='-.', linewidth=0.5, color='grey')

# Save and display the figure
plt.savefig(f"{atoms_number}/energy_conservation_vs_cores.png", dpi=300)
plt.show()

# Clear the current figure for the next plot
plt.clf()

# Second Plot: Total Energy vs Time Steps Across Various Core Counts
energy_columns = ["iteration", "total_energy", "potential_energy"]
core8_data = pd.read_csv(f"{atoms_number}/8_energies.csv", usecols=energy_columns)
core4_data = pd.read_csv(f"{atoms_number}/4_energies.csv", usecols=energy_columns)
core2_data = pd.read_csv(f"{atoms_number}/2_energies.csv", usecols=energy_columns)
core1_data = pd.read_csv(f"{atoms_number}/1_energies.csv", usecols=energy_columns)

# Plotting total energy for various core counts with distinct lines and colors
plt.plot(core8_data['iteration'], core8_data['total_energy'], linestyle='-', color='midnightblue', linewidth=1.8)
plt.plot(core4_data['iteration'], core4_data['total_energy'], linestyle='--', color='gold', linewidth=1.8)
plt.plot(core2_data['iteration'], core2_data['total_energy'], linestyle='-.', color='forestgreen', linewidth=1.8)
plt.plot(core1_data['iteration'], core1_data['total_energy'], linestyle=':', color='tomato', linewidth=1.8)

# Customize labels, title, and legend
plt.legend(["8 Cores", "4 Cores", "2 Cores", "1 Core"], loc='upper right', fontsize=11)
plt.xlabel("Time Steps (fs)", fontsize=13, fontweight='medium', color='black')
plt.ylabel("Total Energy (eV)", fontsize=13, fontweight='medium', color='black')
plt.title("Total Energy Across Core Counts Over Time", fontsize=15, fontweight='bold', color='darkblue')

# Adding a custom grid with different line style
plt.grid(True, linestyle='--', linewidth=0.6, color='lightgray')

# Save the figure and show it
plt.savefig(f"{atoms_number}/total_energy_vs_time_steps_custom.png", dpi=300)
plt.show()

# Clear the figure
plt.clf()