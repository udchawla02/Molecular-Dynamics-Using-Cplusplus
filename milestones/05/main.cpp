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
    try {
        // Initialize atoms parameters
        double mass = 1.0;
        double sigma = 1.0;
        double epsilon = 1.0;
        double time = sqrt(mass * pow(sigma, 2) / epsilon);
        double total_time = 100 * time; // Further reduced total time for debugging
        std::cout << "total_time: " << total_time << std::endl;
        double time_step = time / 1000.0;
        std::cout << "time_step: " << time_step << std::endl;
        double relaxation_time_multiplier = 10.0;
        double relaxation_time;     // relaxation time
        double desired_temperature; // desired temperature

        // Output directory and file paths
        fs::path output_dir = fs::current_path() / "output" / "milestone_05";
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
        std::vector<int> atom_counts = {27, 64, 125}; // Reduced sizes for debugging

        // Vary the number of atoms and record simulation time
        for (int atoms_number : atom_counts) {
            // Calculate dimensions of the cubic lattice
            int n = static_cast<int>(std::cbrt(atoms_number));
            if (n * n * n != atoms_number) {
                throw std::runtime_error("Atom count " + std::to_string(atoms_number) + " is not a perfect cube.");
            }

            // Initialize atoms on a cubic lattice
            double lattice_constant = sigma * 0.8;
            Atoms atoms = lattice(n, n, n, lattice_constant);
            std::cout << "atoms.nb_atoms(): " << atoms.nb_atoms() << std::endl;
            Energy energy(atoms, epsilon, sigma, mass); // initialize energy class

            auto start = std::chrono::high_resolution_clock::now();

            // Main simulation loop
            for (int i = 0; i < static_cast<int>(total_time / time_step); ++i) {
                if (i % 10 == 0) {
                    write_xyz(traj, atoms); // write trajectory
                }
                if (i % 100 == 0) {
                    std::cout << "Step: " << i << " for " << atoms_number << " atoms" << std::endl;
                }
                double old_total_energy = energy.get_total_energy(); // store old total energy
                verlet_step1(atoms, time_step, mass); // update positions
                energy.energy_update(atoms, epsilon, sigma); // update energies
                verlet_step2(atoms, time_step, mass); // update velocities

                // Thermal bathing
                desired_temperature = (i == 0) ? energy.get_temperature() : desired_temperature; // set desired temperature to initial temperature
                if (abs(energy.get_total_energy() - old_total_energy) <= 0.01) { // reached equilibrium point, decrease coupling constant
                    relaxation_time_multiplier = 50;
                }
                relaxation_time = relaxation_time_multiplier * time_step; // relaxation time
                energy.berendsen_thermostat(atoms, desired_temperature, time_step, relaxation_time); // apply Berendsen thermostat
            }

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> simulation_time = end - start;
            std::cout << "Simulation time for " << atoms_number << " atoms: " << simulation_time.count() << " seconds" << std::endl;

            // Log simulation time vs atoms number to CSV
            csv_file << atoms_number << "," << simulation_time.count() << std::endl;
        }

        traj.close();
        csv_file.close();
        std::cout << "Simulation finished successfully!!" << std::endl;

    } catch (const std::runtime_error& e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }

    return 0;
}
