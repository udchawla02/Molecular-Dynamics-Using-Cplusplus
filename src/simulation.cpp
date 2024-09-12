#include "../../../include/Molecular-Dynamics/simulation.h"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cmath> // for std::isnan

Simulation::Simulation() : Energy() {}
Simulation::Simulation(Atoms &new_atoms) : atoms_{new_atoms}, Energy(atoms_, epsilon, sigma, mass) {
    energy_update(atoms_, epsilon, sigma); // Update energy class
    neighbor_list_ = NeighborList(cutoff_radius); // NeighborList object for cutoff
    std::cout << "Initialized Simulation parameters" << std::endl;
}
Simulation::~Simulation() {}

// Helper function to create directory
void create_directory(const std::string &path) {
    std::error_code ec;
    std::filesystem::create_directories(path, ec);
    if (ec) {
        std::cerr << "Error creating directory: " << ec.message() << std::endl;
    } else {
        std::cout << "Directory created successfully: " << path << std::endl;
    }
}

// Helper function to export energy data
void export_energy_data(int step, const std::string &directory, double total_energy, double potential_energy) {
    std::string energy_file_path = directory + "/initial_energy_" + std::to_string(step) + ".csv";
    std::ofstream energy_file(energy_file_path);
    if (energy_file.is_open()) {
        energy_file << "step,total_energy,potential_energy\n";
        energy_file << step << "," << total_energy << "," << potential_energy << "\n";
        energy_file.close();
        std::cout << "Debug: Exporting energy data at step " << step << " to " << energy_file_path << std::endl;
    } else {
        std::cerr << "Error opening energy file: " << energy_file_path << std::endl;
    }
}

// Helper function to export melting point data
void export_melting_point_data(const std::string &base_dir, int cluster_size, double melting_point, int iterations, double added_energy, int atoms_number, double latent_energy) {
    std::string melting_point_file_path = base_dir + "/melting_points_vs_cluster_size.csv";
    bool file_exists = std::filesystem::exists(melting_point_file_path);
    std::ofstream melting_point_file(melting_point_file_path, std::ios_base::app); // Append mode

    if (melting_point_file.is_open()) {
        if (!file_exists) {
            melting_point_file << "cluster_size,melting_point,number_of_iterations_till_melting,added_energy,atoms_number,latent_energy\n";
        }
        melting_point_file << cluster_size << "," << melting_point << "," << iterations << "," << added_energy << "," << atoms_number << "," << latent_energy << "\n";
        melting_point_file.close();
        std::cout << "Debug: Exporting melting point data for cluster size " << cluster_size << " to " << melting_point_file_path << std::endl;
    } else {
        std::cerr << "Error opening melting point file: " << melting_point_file_path << std::endl;
    }
}

// Main simulation loop
void Simulation::initial_loop() {
    bool equilibrium = false;
    double previous_energy = get_total_energy();
    double melting_point = 0.0;
    bool melting_detected = false;
    int melting_iteration = 0;

    for (int i = 0; i < total_steps; ++i) {
        std::cout << "Step: " << i << std::endl;
        if (equilibrium) std::cout << "Equilibrium reached" << std::endl;
        std::cout << "Initial steps: " << i << "  Current temp: " << get_temperature() << "  Current potential: " << get_potential_energy() << "  Total energy: " << get_total_energy() << std::endl;

        double current_energy = get_total_energy();

        // Detect melting point based on significant energy change
        if (!melting_detected && std::abs(current_energy - previous_energy) > 0.1) {
            melting_point = get_temperature();
            melting_iteration = i;
            melting_detected = true;
        }

        previous_energy = current_energy;

        verlet_step1(atoms_, time_step, mass);
        neighbor_list_.update(atoms_);
        update_gupta(atoms_, neighbor_list_, cutoff_radius);
        verlet_step2(atoms_, time_step, mass);

        // Thermal bathing
        if (i == stop_thermostate_after_steps) {
            equilibrium = true;
            relaxation_time_multiplier = relaxation_time_multiplier_final_value; // This should be big enough to reduce thermostat effect
        }
        relaxation_time = relaxation_time_multiplier * time_step; // Relaxation time
        berendsen_thermostat(atoms_, desired_temperature, time_step, relaxation_time);

        // Check for NaN values
        if (std::isnan(get_total_energy()) || std::isnan(get_potential_energy()) || std::isnan(get_temperature())) {
            std::cerr << "Error: NaN detected at step " << i << std::endl;
            MPI_Finalize();
            exit(1);
        }

        if (i % 10 == 0) {
            std::string step_directory = directory + "/step_" + std::to_string(i);
            create_directory(step_directory);

            std::cout << "Debug: Exporting XYZ at step " << i << std::endl;
            export_xyz_initial(step_directory, i, atoms_); // Export to XYZ file
            export_energy_data(i, step_directory, get_total_energy(), get_potential_energy());
        }
    }

    if (melting_detected) {
        int atoms_number = static_cast<int>(atoms_.positions.rows());
        double latent_energy = get_total_energy() - previous_energy;
        export_melting_point_data(directory, layer_numbers, melting_point, melting_iteration, add_energy, atoms_number, latent_energy);
    }
}

// Overloaded initial loop with Domain
void Simulation::initial_loop(Domain domain) {
    bool equilibrium = false;
    double previous_energy = get_total_energy();
    double melting_point = 0.0;
    bool melting_detected = false;
    int melting_iteration = 0;

    for (int i = 0; i < total_steps; ++i) {
        std::cout << "Step: " << i << std::endl;
        double old_total_energy = get_total_energy();
        verlet_step1(atoms_, time_step, mass);
        domain.exchange_atoms(atoms_); // Exchange atoms between domains after updating positions
        domain.update_ghosts(atoms_, cutoff_radius * 2); // Update ghost atoms before calculating forces
        neighbor_list_.update(atoms_);
        update_gupta(atoms_, neighbor_list_, cutoff_radius, domain);
        verlet_step2(atoms_, time_step, mass);

        // Detect melting point based on significant energy change
        double current_energy = get_total_energy();
        if (!melting_detected && std::abs(current_energy - previous_energy) > 0.1) {
            melting_point = get_temperature();
            melting_iteration = i;
            melting_detected = true;
        }
        previous_energy = current_energy;

        // Thermal bathing
        if (i == stop_thermostate_after_steps) {
            relaxation_time_multiplier = relaxation_time_multiplier_final_value; // This should be big enough to reduce thermostat effect
        }
        if (std::abs(get_total_energy() - old_total_energy) < 0.01) {
            domain.disable(atoms_);
            std::cout << "Equilibrium reached, exiting now!" << std::endl;
            MPI_Finalize();
            exit(1);
        }
        relaxation_time = relaxation_time_multiplier * time_step; // Relaxation time
        berendsen_thermostat(atoms_, desired_temperature, time_step, relaxation_time);

        // Check for NaN values
        if (std::isnan(get_total_energy()) || std::isnan(get_potential_energy()) || std::isnan(get_temperature())) {
            std::cerr << "Error: NaN detected at step " << i << std::endl;
            MPI_Finalize();
            exit(1);
        }

        if (i % 10 == 0) {
            std::string step_directory = directory + "/step_" + std::to_string(i);
            create_directory(step_directory);

            std::cout << "Debug: Exporting XYZ at step " << i << std::endl;
            domain.disable(atoms_);
            export_xyz_initial(step_directory, i, atoms_);
            domain.enable(atoms_);
            domain.exchange_atoms(atoms_);
            domain.update_ghosts(atoms_, cutoff_radius * 2);
            export_energy_data(i, step_directory, get_total_energy(), get_potential_energy());
        }
    }

    if (melting_detected) {
        int atoms_number = static_cast<int>(atoms_.positions.rows());
        double latent_energy = get_total_energy() - previous_energy;
        export_melting_point_data(directory, layer_numbers, melting_point, melting_iteration, add_energy, atoms_number, latent_energy);
    }
}

void Simulation::relaxation_loop(int iteration) {
    std::cout << "--------------------Starting relaxation loop number: " << iteration << "---------------------------" << std::endl;
    double total_temp = 0.0;
    for (int i = 0; i < relaxation_steps; ++i) {
        std::cout << "Relaxation steps: " << i << "  Current temp: " << get_temperature() << "  Current potential: " << get_potential_energy() << "  Total energy: " << get_total_energy() << std::endl;
        total_temp += get_temperature();
        verlet_step1(atoms_, time_step, mass);
        neighbor_list_.update(atoms_);
        update_gupta(atoms_, neighbor_list_, cutoff_radius);
        verlet_step2(atoms_, time_step, mass);

        // Check for NaN values
        if (std::isnan(get_total_energy()) || std::isnan(get_potential_energy()) || std::isnan(get_temperature())) {
            std::cerr << "Error: NaN detected at step " << i << std::endl;
            MPI_Finalize();
            exit(1);
        }

        if (i % 10 == 0) {
            std::string step_directory = directory + "/relax_" + std::to_string(iteration) + "_step_" + std::to_string(i);
            create_directory(step_directory);
            export_xyz_relax(step_directory, iteration * relaxation_steps + i, atoms_);
        }
    }
    double average_temperature = total_temp / relaxation_steps;
    export_data(iteration, energy_file, get_total_energy(), average_temperature, get_potential_energy(), continue_old_experiment);
}

void Simulation::add_heat() {
    deposit_heat(atoms_, add_energy);
}