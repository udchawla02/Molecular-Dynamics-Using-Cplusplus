#include "../../../include/Molecular-Dynamics/verlet.h"
#include "../../../include/Molecular-Dynamics/xyz.h"
#include "../../../include/Molecular-Dynamics/energy.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <cmath>

namespace fs = std::filesystem;

int main() {
    try {
        // Initialize atoms parameters
        double mass = 1.0;
        double sigma = 1.0;
        double epsilon = 1.0;
        double time = sqrt(mass * pow(sigma, 2) / epsilon);
        double total_time = 5000 * time;

        // Construct the path to the lj54.xyz file
        fs::path file_path = fs::current_path() / "milestones" / "04" / "lj54.xyz";
        std::string filename = file_path.string();

        std::cout << "Attempting to open file: " << filename << std::endl;

        // Check if the file exists
        if (!fs::exists(file_path)) {
            throw std::runtime_error("Error: File does not exist: " + filename);
        }

        std::ifstream file_check(filename);
        if (!file_check.is_open()) {
            throw std::runtime_error("Error: Could not open file " + filename);
        }
        file_check.close();

        // Reading initial positions and velocities from xyz file
        auto [names, initial_positions, initial_velocities] = read_xyz_with_velocities(filename);

        // Creating visualization file milestones/04/output
        fs::path output_dir = fs::current_path() / "output" / "milestone_04";
        if (!fs::exists(output_dir)) {
            fs::create_directories(output_dir);
            std::cout << "Created directory: " << output_dir << std::endl;
        }

        std::ofstream traj(output_dir / "traj.xyz");
        if (!traj.is_open()) {
            throw std::runtime_error("Error: Could not create traj.xyz output file");
        }

        std::ofstream energy_file(output_dir / "energies_vs_time_step.csv");
        if (!energy_file.is_open()) {
            throw std::runtime_error("Error: Could not create energies_vs_time_step.csv output file");
        }
        energy_file << "time_step,total_energy,kinetic_energy,potential_energy" << std::endl;

        // Main simulation loop
        for (double time_step = 1e-4; time_step < 30e-3; time_step += 1e-3) {
            // initializing atoms with positions and initial velocities
            Atoms atoms{initial_positions, initial_velocities};
            Energy energy(atoms, epsilon, sigma, mass); // initialize energy class
            std::cout << "time_step: " << time_step << std::endl;
            for (int i = 0; i < total_time; ++i) {
                if (i % 10 == 0) write_xyz(traj, atoms);    // writing to xyz file
                double old_total_energy = energy.get_total_energy();    // storing the old total energy
                verlet_step1(atoms, time_step, mass);   // updating positions
                energy.energy_update(atoms, epsilon, sigma);    // updating energies
                verlet_step2(atoms, time_step, mass);   // updating velocities
                // check if the system reached equilibrium point to move on to the next time step
                if (abs(energy.get_total_energy() - old_total_energy) < 0.001 && i > 2) {
                    std::cout << "Equilibrium reached for this step, exiting now!!" << std::endl;
                    energy_file << time_step << "," << energy.get_total_energy() << "," << energy.get_kinetic_energy()
                                << "," << energy.get_potential_energy() << std::endl;
                    break;
                }
            }
        }
        traj.close();
        energy_file.close();
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
