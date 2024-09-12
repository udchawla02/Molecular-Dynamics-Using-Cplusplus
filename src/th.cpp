#include "../../../include/Molecular-Dynamics/th.h"

double temperature_cur(const Atoms atoms) {
    // Calculate the kinetic energy and return temperature
    double kinetic_energy = 0.0;

    // Loop over all atoms and calculate kinetic energy manually
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        Eigen::Vector3d velocity = atoms.velocities.col(i);  // Get the velocity vector
        kinetic_energy += 0.5 * atoms.masses(i) * velocity.dot(velocity);  // Use dot product to calculate squared norm
    }

    // Return the temperature
    return kinetic_energy / (KB_EV_K_3_2 * atoms.nb_atoms());
}

void berendsen_thermostat(Atoms &atoms, double temperature, double dt, double relaxation_time) {
    double temp_cur = temperature_cur(atoms);
    double lambda = temp_cur > 0. ? sqrt(1 + (temperature / temp_cur - 1) * dt / relaxation_time) : 1.0;
    atoms.velocities *= lambda;
}

void berendsen_thermostat_decomposed(Atoms &atoms, double temperature, double dt, double relaxation_time, size_t nb_local, double mass) {
    double kinetic_energy = 0.0;

    // Compute kinetic energy of the first `nb_local` atoms manually
    for (int i = 0; i < nb_local; ++i) {
        Eigen::Vector3d velocity = atoms.velocities.col(i);  // Get the velocity vector
        kinetic_energy += 0.5 * mass * velocity.dot(velocity);  // Use dot product to calculate squared norm
    }

    double temp_cur = kinetic_energy / (KB_EV_K_3_2 * nb_local);
    double lambda = temp_cur > 0. ? sqrt(1 + (temperature / temp_cur - 1) * dt / relaxation_time) : 1.0;

    atoms.velocities *= lambda;
}
