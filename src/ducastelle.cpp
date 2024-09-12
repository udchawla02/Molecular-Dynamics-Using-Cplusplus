/*
 * Copyright 2021 Lars Pastewka
 *
 * ### MIT license
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include "../include/Molecular-Dynamics/ducastelle.h"

// Existing ducastelle function with Atoms& as the first argument
double ducastelle(Atoms &atoms, const NeighborList &neighbor_list, double cutoff, double A, double xi, double p, double q, double re) {
    auto cutoff_sq{cutoff * cutoff};
    double xi_sq{xi * xi};

    atoms.forces.setZero();  // Reset energies and forces

    // Compute embedding energies
    Eigen::ArrayXd embedding(atoms.nb_atoms());
    embedding.setZero();

    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            Eigen::Vector3d distance_vector{atoms.positions.col(i) - atoms.positions.col(j)};
            auto distance_sq = distance_vector.squaredNorm();
            if (distance_sq < cutoff_sq) {
                double density_contribution{xi_sq * std::exp(-2 * q * (std::sqrt(distance_sq) / re - 1.0))};
                embedding(i) += density_contribution;
                embedding(j) += density_contribution;
            }
        }
    }

    embedding = -embedding.sqrt();
    Eigen::ArrayXd energies{embedding};

    // Compute forces
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            double d_embedding_density_i = embedding(i) != 0 ? 1 / (2 * embedding(i)) : 0;
            Eigen::Vector3d distance_vector{atoms.positions.col(i) - atoms.positions.col(j)};
            auto distance_sq = distance_vector.squaredNorm();
            if (distance_sq < cutoff_sq) {
                double distance = std::sqrt(distance_sq);
                double d_embedding_density_j = embedding(j) != 0 ? 1 / (2 * embedding(j)) : 0;

                double repulsive_energy = 2 * A * std::exp(-p * (distance / re - 1.0));
                double d_repulsive_energy = -repulsive_energy * p / re;

                double fac = -2 * q / re * xi_sq * std::exp(-2 * q * (distance / re - 1.0));

                Eigen::Array3d pair_force = (d_repulsive_energy + fac * (d_embedding_density_i + d_embedding_density_j)) * distance_vector.normalized();

                repulsive_energy *= 0.5;
                energies(i) += repulsive_energy;
                energies(j) += repulsive_energy;

                atoms.forces.col(i) -= pair_force;
                atoms.forces.col(j) += pair_force;
            }
        }
    }

    return energies.sum();
}

// Overloaded version of ducastelle that accepts an int as the first argument
double ducastelle(int nb_local, Atoms &atoms, const NeighborList &neighbor_list, double cutoff, double A, double xi, double p, double q, double re) {
    double total_energy = 0.0;

    // Loop over local atoms and apply the potential calculation (example)
    for (int i = 0; i < nb_local; ++i) {
        // Here you would apply the same logic as in the existing ducastelle function
        // for atoms[i] and their interactions within the neighbor list.
        total_energy += ducastelle(atoms, neighbor_list, cutoff, A, xi, p, q, re);
    }

    return total_energy;
}
