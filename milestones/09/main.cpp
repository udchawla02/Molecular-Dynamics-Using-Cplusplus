#include "../../../include/Molecular-Dynamics/atoms.h"
#include "mpi.h"
#include "../../../include/Molecular-Dynamics/mpi_support.h"
#include "../../../include/Molecular-Dynamics/verlet.h"
#include "../../../include/Molecular-Dynamics/xyz.h"
#include <../../../include/Molecular-Dynamics/domain.h>
#include <../../../include/Molecular-Dynamics/ducastelle.h>
#include <iomanip>
#include <iostream>
#include <../../../include/Molecular-Dynamics/neighbors.h>
#include <../../../include/Molecular-Dynamics/th.h>

// simulation settings
const double DT{5};
const double TAU{1000.};
const size_t SCALE_AFTER_N{(size_t)round(10 * TAU)};
const double TEMPERATURE{200.};
const double STRAIN_RATE_HZ{1e8};
const double STRETCH_BY_PERCENT{30};
const std::string FILENAME{"whisker_large_cold"};

// other constants
const double CUTOFF_RADIUS{10.};
const size_t TIMESTEPS{(size_t)ceil((STRETCH_BY_PERCENT / 100.) /
                                    (STRAIN_RATE_HZ) / (DT * 1e-15))};
const size_t PLOT_EVERY{TIMESTEPS / 1000};
const double SPACING{4.079 / sqrt(2)};
const double BORDER{50 * SPACING};

// Vec3_t should be Eigen::Vector3d
using Vec3_t = Eigen::Vector3d;

/// @brief Subtract the mean velocity of the atoms in a system from it to
/// prevent bulk movement without affecting relative velocities, making
/// temperature measurements more accurate.
/// @param vel the velocities to correct
void subtract_mean(Velocities_t &vel) {
    Vec3_t v_avg{vel.rowwise().mean()};
    vel.row(0) -= v_avg(0);
    vel.row(1) -= v_avg(1);
    vel.row(2) -= v_avg(2);
}

/// @brief Calculate the volume of the atoms provided and move them to fit
/// neatly in the volume. Leave a lot of extra space in the non-periodic x- and
/// y-direciton and just a little bit in the z-direction.
/// @param atoms the atoms to move and construct a bounding volume of
/// @return the volume bounding the atoms as an (x,y,z)-vector
Vec3_t calculate_volume(Atoms &atoms) {
    // move the atoms such that the lowest z-value is zero and the volume
    // specified later surrounds the whisker in x,y directions
    Vec3_t min_pos{atoms.positions.rowwise().minCoeff()};
    atoms.positions.row(0) -= min_pos(0);
    atoms.positions.row(1) -= min_pos(1);
    atoms.positions.row(2) -= min_pos(2);
    min_pos = atoms.positions.rowwise().minCoeff();
    Vec3_t max_pos{atoms.positions.rowwise().maxCoeff()};
    Vec3_t volume{max_pos - min_pos};
    atoms.positions.row(0) += BORDER;
    atoms.positions.row(1) += BORDER;
    // increase the volume in x and y direction as much as desired so atoms
    // don't go out of bounds, only increase in z-direction by some small
    // epsilon so atoms are not lost when the domain is activated.
    volume(0) += 2.0 * BORDER;
    volume(1) += 2.0 * BORDER;
    volume(2) += SPACING / 4.0;
    return volume;
}

/// @brief Print information about the current state of the system to the
/// specified output stream in csv format
/// @param out the output stream, may be `std::cout` or an `ofstream`
/// @param atoms the current state of the system
/// @param volume the current volume of the simulation domain
/// @param force_zz_g the global force in zz direction
/// @param e_pot_global the global potential energy
/// @param t the current time in femtoseconds
/// @param v_z_0 the initial volume in z-direction to divide strain by
void print_results(std::ostream &out, const Atoms &atoms, const Vec3_t volume,
                   const double force_zz_g, const double e_pot_global,
                   const double t, const double v_z_0) {
    // header:
    // "temp,f_zz,v_z,e_pot,t,v_z_0,rate,v_x,v_y"
    out << std::fixed << std::setprecision(15)

        << temperature_cur(atoms) << "," << force_zz_g << "," << volume(2)
        << "," << e_pot_global << "," << t << "," << v_z_0 << ","
        << STRAIN_RATE_HZ << "," << volume(0) << "," << volume(1)

        << std::endl;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // initialize gold cluster from file
    auto [names, positions,
          velocities]{read_xyz_with_velocities(FILENAME + ".xyz")};
    Atoms atoms{Atoms(names, positions)};
    // miniscule random velocities to get kinetic energy in all modes
    atoms.velocities.setRandom();
    atoms.velocities *= 1e-4;
    // deduct mean velocity to prevent bulk movement
    subtract_mean(atoms.velocities);

    // open files
    // Commenting this out since `ELEM_NAME_TO_MASS` is undefined
    // double m{ELEM_NAME_TO_MASS.at("Au")};
    double m = 197.0; // Placeholder mass for gold (you can use this temporarily)

    std::ofstream traj(FILENAME + "_traj.xyz");
    std::ofstream results(FILENAME + ".csv");

    // create neighbour list and initialize forces
    NeighborList neighbour_list;
    (void)ducastelle(atoms, neighbour_list);

    // find the number of MPI processes
    int nb_processes{0};
    MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);

    // calculate the volume of the domain
    Vec3_t volume{calculate_volume(atoms)};

    // here, only decompose along the z-axis,  where the whisker extends and
    // use periodic boundary in z-direction
    Domain domain(MPI_COMM_WORLD, volume, {1, 1, nb_processes}, {0, 0, 1});

    if (domain.rank() == 0) {
        // write csv header to stdout
        std::cout << "temp,f_zz,v_z,e_pot,t,v_z_0,rate,v_x,v_y" << std::endl;
        results << "temp,f_zz,v_z,e_pot,t,v_z_0,rate,v_x,v_y" << std::endl;
    }

    // switch to decomposed state
    domain.enable(atoms);

    // initial force calculation
    domain.update_ghosts(atoms,
                         CUTOFF_RADIUS *
                             2.); // border radius is twice the cutoff radius
    (void)ducastelle(atoms, neighbour_list);

    // main simulation loop
    double t{0};
    double force_zz_local{0.};
    double v_z_0{volume(2)};
    for (size_t i{0}; i < TIMESTEPS + SCALE_AFTER_N; i++) {
        if (i >= SCALE_AFTER_N) {
            if (i == SCALE_AFTER_N) {
                domain.disable(atoms);
                volume = calculate_volume(atoms);
                v_z_0 = volume(2);
                domain.enable(atoms);
            }
            // continuously rescale the domain to stretch it in z-direction
            volume(2) = v_z_0 + v_z_0 * STRAIN_RATE_HZ * t * (1e-15);
            domain.scale(atoms, volume);
        }

        if (i >= SCALE_AFTER_N && (i - SCALE_AFTER_N) % PLOT_EVERY == 0) {
            double force_zz_g{MPI::allreduce(force_zz_local / PLOT_EVERY,
                                             MPI_SUM, domain.communicator())};
            force_zz_local = 0.0;

            domain.disable(atoms);
            subtract_mean(atoms.velocities);

            if (domain.rank() == 0) {
                write_xyz(traj, atoms);
                std::ofstream latest("latest_whisker.xyz");
                write_xyz(latest, atoms);
                latest.close();
            }
            domain.enable(atoms);
            domain.update_ghosts(atoms, CUTOFF_RADIUS * 2.);
            (void)ducastelle(atoms, neighbour_list);
        }

        if (i >= SCALE_AFTER_N) {
            t += DT;
        }
    }

    if (domain.rank() == 0)
        std::cout << std::endl;
    MPI_Finalize();
    return 0;
}
