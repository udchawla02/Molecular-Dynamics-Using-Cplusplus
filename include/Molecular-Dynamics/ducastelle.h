#ifndef YAMD_DUCASTELLE_H
#define YAMD_DUCASTELLE_H

#include "atoms.h"
#include "neighbors.h"

// Existing function declaration
double ducastelle(Atoms &atoms, const NeighborList &neighbor_list, double cutoff = 10.0, double A = 0.2061,
                  double xi = 1.790, double p = 10.229, double q = 4.036, double re = 4.079 / sqrt(2));

// Overloaded function to accept an integer (nb_local)
double ducastelle(int nb_local, Atoms &atoms, const NeighborList &neighbor_list, double cutoff = 10.0, double A = 0.2061,
                  double xi = 1.790, double p = 10.229, double q = 4.036, double re = 4.079 / sqrt(2));

#endif //YAMD_GUPTA_H
