#ifndef EXPORT_H
#define EXPORT_H

#include <string>
#include <fstream>

void export_data(const std::string& filename, double total_energy, double potential_energy);
void export_data(std::ofstream& energy_file, double total_energy, double average_temperature, double potential_energy, bool continue_old_experiment);

#endif // EXPORT_H
