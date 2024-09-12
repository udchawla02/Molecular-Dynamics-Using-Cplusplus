/*
* Copyright 2021 Robert Schütze
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
* 
* (Created by Robert Schütze on 27.08.2021.)
*/

#include <gtest/gtest.h>
#include "../include/Molecular-Dynamics/verlet.h"
#include "../include/Molecular-Dynamics/energy.h"
#include "../include/Molecular-Dynamics/lj_direct_summation.h"

// launch single atom with v = (1, 1, 1), bring it to target temperature = 2.0
// and check whether it follows the desired exponential curve
TEST(BerendsenThermostatTest, SingleAtomExponentialDecay) {
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;

    double dt = 0.01; // timestep
    double tau = 1000.*dt; // big enough relaxation time for the approximation e^-dt/tau=-dt/tau to hold
    int steps = 10000; // number of simulation steps
    double T_target = 0.1; // target temperature

    Atoms atoms(lattice(1, 1, 1, 1));
    atoms.velocities.setOnes();
    Energy energy(atoms, epsilon, sigma, mass);

    double T_init, T_curr, T_calc;

    for(int i=1; i<=steps; i++){ // main sim loop
        verlet_step1(atoms, dt,1);
        energy.energy_update(atoms, epsilon, sigma);
        verlet_step2(atoms, dt,1);
        T_init = energy.get_temperature();
        energy.berendsen_thermostat(atoms, T_target, dt, tau);
        T_curr = energy.get_temperature();
        T_calc = T_target+(T_init-T_target)*exp(-dt/tau);
        EXPECT_NEAR(T_curr, T_calc, T_calc*0.01); // check whether T_curr=T_calc+-1%
    }
}

// equilibrate an atom cluster to T = 0.1 epsilon
// and check whether T converges to the desired T
TEST(BerendsenThermostatTest, FinalConvergence) {
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;

    double dt = 0.01; // timestep
    double tau = 10.*dt; // relaxation time
    int steps = 10000; // number of simulation steps
    double T_target = 1.0;

    Atoms atoms(lattice(5, 5, 5, sigma*1.12)); // initialize 5x5x5 lattice
    Energy energy(atoms, epsilon, sigma, mass);
    for(int i=1; i<=steps; i++){ // main sim loop
        std::cout << "step " << i << std::endl; // print current step
        verlet_step1(atoms, dt, 1);
        energy.energy_update(atoms, epsilon, sigma);
        verlet_step2(atoms, dt, 1);
        energy.berendsen_thermostat(atoms, T_target, dt, tau);
        if(i>steps/2){ // only check after convergence
//            std::cout << energy.get_temperature() << std::endl; // print current temperature
            EXPECT_NEAR(energy.get_temperature(), T_target, T_target*0.2); // check whether T=T_target+-20%
        }
    }
}