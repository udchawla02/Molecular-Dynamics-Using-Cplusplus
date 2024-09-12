#include "../../../include/Molecular-Dynamics/simulation_data.h"

SimulationData::SimulationData() {
    // Initialize all simulation parameters
//    cluster_name = "cluster_923"; // layer_number = 6
    layer_numbers = 7;
    atomic_distance = 2.885; // atomic distance from reference clusters - corresponds to 408 pm lattice constant
    mass = 196.9665* 103.6; // atomic mass of Gold (https://www.nuclear-power.com/gold-atomic-number-mass-density/)
    total_steps = 5;
    time_step = 1.0; // time step in fs
    cutoff_radius = 10.0;    // cutoff radius for EAM potential
    relaxation_time_multiplier = 10; // relaxation time = relaxation time multiplier * time_step in fs
    stop_thermostate_after_steps = 500; // stop thermostat after this number of steps
    relaxation_time_multiplier_final_value = 1e10; // after the system arrives at desired temp (stop thermostat and relax system)
    desired_temperature = 500.0; // desired temperature (only in the start) in K

    // Relaxation experiment parameters
    relaxation_steps = 2; // number of relaxation steps
    expermint_num = 2;    // number of experiments
    add_energy = 0.01;    // energy added in each experiment

    // MPI parameters
    domain_length = {30.0, 30.0, 30.0}; // domain length in Angstrom
    domain_grid = {2, 2, 2}; // number of domains in each direction
    domain_periodicity = {0, 0, 0}; // periodicity of the domain in each direction

    // choose whether to continue old experiment or not
    continue_old_experiment = false;
    old_experiment_file = "traj_1415_4990_initial.xyz"; // name of the old experiment file (ONLY if continuing old experiment)

    create_directories_and_files();
}

void SimulationData::create_directories_and_files() {
    const char* base_dir_env = std::getenv("BASE_DIR");
    std::string base_dir = base_dir_env ? base_dir_env : ".";
    std::string number_of_layers = std::to_string(layer_numbers);

    std::string milestone_number = "08";
    directory = base_dir + "/output/milestone_" + milestone_number + "/" + number_of_layers + "/";

    std::error_code ec;
    std::filesystem::create_directories(directory, ec);
    if (ec) {
        std::cerr << "Error creating directory: " << ec.message() << std::endl;
    } else {
        std::cout << "Directory created successfully: " << directory << std::endl;
    }

    energy_file = std::ofstream(directory + number_of_layers + "_energies.csv");
    if (!energy_file) {
        std::cerr << "Error opening energy file: " << directory + number_of_layers + "_energies.csv" << std::endl;
    } else {
        std::cout << "Energy file created successfully: " << directory + number_of_layers + "_energies.csv" << std::endl;
    }

    old_experiment_file = directory + old_experiment_file;
}