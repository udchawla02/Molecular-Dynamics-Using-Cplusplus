#include "../include/Molecular-Dynamics/export.h"
#include <fstream>

void export_data(const std::string& filename, double total_energy, double potential_energy) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "Total Energy,Potential Energy\n";
        file << total_energy << "," << potential_energy << "\n";
        file.close();
    }
}

void export_data(std::ofstream& energy_file, double total_energy, double average_temperature, double potential_energy, bool continue_old_experiment) {
    if (energy_file.is_open()) {
        energy_file << "Iteration,Total Energy,Average Temperature,Potential Energy,Continue Old Experiment\n";
        energy_file << total_energy << "," << average_temperature << "," << potential_energy << "," << continue_old_experiment << "\n";
    }
}
