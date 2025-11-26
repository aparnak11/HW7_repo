#include <iostream>
#include <fstream>          // for file writing
//#include "chamber_heat.h"

bool is_thruster_on(int k, double dt)
{
    // pulse parameters
    double pulse_start_time = 0.0; // pulse start time in seconds
    double pulse_duration = 0.15; // pulse end time in seconds
    double pulse_period = 0.2; // pulse period in seconds

    bool thruster_on = false;
    double current_time = k * dt;
    double time_in_cycle = fmod(current_time, pulse_period);
    if (time_in_cycle >= pulse_start_time && time_in_cycle < (pulse_start_time + pulse_duration))
    {
        thruster_on = true; // thruster is on during this time step
    }
    return thruster_on;
}

int main() {
    // simulation parameters
    double L = 0.01; // length in meters
    double dt = 1e-4; // time step in seconds
    double ncells = 100; // number of cells
    double nsteps = 15000; // number of time steps
    double dx = L / ncells; // cell size in meters

    double A[100][15000];

    for (int i=0; i<ncells; i++) {
        for (int k=0; k<nsteps; k++) {
            A[i][k] = 300.0; // initial temperature in K
        }
    }
    

    // material properties
    double k_cond = 65.0; // thermal conductivity in W/m-K
    double alpha = 1.7e-5; // thermal diffusivity in m^2/s
    double hg = 1.5e4; // convective heat transfer coefficient in W/m^2-K
    double Tog = 2500.0; // hot gas temperature in K


    for (int k=1; k<nsteps; k++) {
        A[0][k] = A[0][k-1] + alpha * dt * (A[1][k] - 2*A[0][k]) / (dx * dx);
        if (is_thruster_on(k, dt)) {
            A[0][k] = -A[0][k-1] + hg/k_cond * (Tog - A[0][k]) * dx;
        }
        printf("Time: %f s, T[0]: %f K\n", k*dt, A[0][k]);
    }


    return 0;
}
