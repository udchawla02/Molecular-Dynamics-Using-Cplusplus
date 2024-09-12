import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = [8.00, 4.50]
plt.rcParams["figure.autolayout"] = False

columns = ["atoms_number", "simulation_time"]
df = pd.read_csv("simulation_time_vs_atoms_number.csv", usecols=columns)

plt.plot(df.atoms_number, df.simulation_time, color='purple', marker='o', linestyle='--', linewidth=2)

plt.xlabel("Atoms Number", fontsize=12, color='darkblue')
plt.ylabel("Simulation Time (s)", fontsize=12, color='darkblue')
plt.title("Simulation Time vs Atoms Number", fontsize=14, color='darkgreen')

plt.grid(True, which='both', linestyle=':', linewidth=0.75, color='gray')

plt.savefig("simulation_time_vs_atoms_number_custom.png")
plt.show()
