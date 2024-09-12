#include "../../../include/Molecular-Dynamics/energy.h"
#include "../../../include/Molecular-Dynamics/verlet.h"
#include "../../../include/Molecular-Dynamics/xyz.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <chrono>
#include <vector>

namespace fs = std::filesystem;

int main() {
    // Initialize atoms parameters
    double mass = 1.0;
    double sigma = 1.0;
    double epsilon = 1.0;
    double time = sqrt(mass * pow(sigma, 2) / epsilon);
    double total_time = 5000 * time;
    double time_step = time / 1000.0;
    double cutoff_radius = 1.5;    // cutoff radius for LJ potential
    double relaxation_time_multiplier = 10.0;
    double relaxation_time;     // relaxation time
    double desired_temperature; // desired temperature

    // Output directory and file paths
    fs::path output_dir = fs::current_path() / "output" / "milestone_06";
    if (!fs::exists(output_dir)) {
        fs::create_directories(output_dir);
        std::cout << "Created directory: " << output_dir << std::endl;
    }

    std::ofstream traj(output_dir / "traj.xyz");
    if (!traj.is_open()) {
        throw std::runtime_error("Error: Could not create traj.xyz output file");
    }

    std::ofstream csv_file(output_dir / "simulation_time_vs_atoms_number.csv");
    if (!csv_file.is_open()) {
        throw std::runtime_error("Error: Could not create simulation_time_vs_atoms_number.csv output file");
    }
    csv_file << "atoms_number,simulation_time" << std::endl;

    // Specific atom counts to test
    std::vector<int> lattice_sizes = {3, 4, 5, 6, 7, 8}; // Adjust sizes as needed

    // Vary the number of atoms and record simulation time
    for (int n : lattice_sizes) {
        // Initialize atoms on a cubic lattice
        double lattice_constant = sigma * 0.8;
        Atoms atoms = lattice(n, n, n, lattice_constant);
        Energy energy(atoms, epsilon, sigma, mass); // initialize energy class
        NeighborList neighbor_list(cutoff_radius); // initialize NeighborList object for cutoff

        auto start = std::chrono::high_resolution_clock::now();

        // Main simulation loop
        for (int i = 0; i < static_cast<int>(total_time); ++i) {
            if (i % 10 == 0) {
                write_xyz(traj, atoms);
            }
            double old_total_energy = energy.get_total_energy();
            verlet_step1(atoms, time_step, mass);
            energy.update_neighbors(atoms, neighbor_list, epsilon, sigma, mass);
            verlet_step2(atoms, time_step, mass);
            // thermal bathing
            // to preserve the temperature, we assume that the temperature is constant, and then we apply the Berendsen thermostat
            desired_temperature = (i == 0) ? energy.get_temperature() : desired_temperature;
            if (abs(energy.get_total_energy() - old_total_energy) <= 0.01) { // reached equilibrium point, decrease coupling constant
                relaxation_time_multiplier = 50;
            }
            relaxation_time = relaxation_time_multiplier * time_step; // relaxation time
            energy.berendsen_thermostat(atoms, desired_temperature, time_step, relaxation_time);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> simulation_time = end - start;
        int atoms_number = n * n * n; // total number of atoms
        std::cout << "Simulation time for " << atoms_number << " atoms: " << simulation_time.count() << " seconds" << std::endl;

        // Log simulation time vs atoms number to CSV
        csv_file << atoms_number << "," << simulation_time.count() << std::endl;
    }

    traj.close();
    csv_file.close();
    std::cout << "Simulation finished successfully!!" << std::endl;

    return 0;
}
