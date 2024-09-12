import matplotlib.pyplot as plt
import numpy as np

# Generating new data for strain and stress force
strain_values = np.linspace(0, 0.30, 100)
stress_force_values = np.random.normal(4.8, 1.2, 100)  # Adjusted values

# Setting up a new figure with a slightly altered size
plt.figure(figsize=(9, 4.5))
plt.plot(strain_values, stress_force_values, linestyle='--', color='blue', linewidth=1.7)

# Customizing the plot with new headings and labels
plt.title("Force vs Strain Curve for Material B", fontsize=15, fontweight='semibold', color='darkblue')
plt.xlabel("Strain (%)", fontsize=13, fontweight='light', color='darkslategray')
plt.ylabel("Stress Force (N)", fontsize=13, fontweight='light', color='darkslategray')

# Altering grid appearance for uniqueness
plt.grid(True, linestyle='-.', color='lightgrey', alpha=0.8)

# Adding distinct annotations with varied descriptive information
plt.text(0.22, 7.5, 'Material ID: B1056, Alloy: Aluminum', fontsize=9, color='brown', style='italic')
plt.text(0.22, 7.1, 'Temp: 310K', fontsize=9, color='brown', style='italic')
plt.text(0.22, 6.7, 'Trial Number: 12', fontsize=9, color='brown', style='italic')

# Save the plot with a distinct name and display it
plt.savefig('custom_stress_vs_strain_curve.png', dpi=300)
plt.show()