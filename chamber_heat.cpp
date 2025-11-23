#include <iostream>

int main() {
    double domain_length = 0.01; // Length of the domain in meters
    double dt = 1e-4; // Time step in seconds
    int num_cells = 100; // Number of spatial cells
    int num_time_steps = 15000; // Number of time steps

    // material properties
    double k = 65.0; // Thermal conductivity in W/(m·K)
    double alpha = 1.7e-5; // Thermal diffusivity in m^2/s
    double hg = 1.5e4; // Convective heat transfer coefficient in W/(m^2·K)
    double T0g = 2500; // Gas stagnation temperature in Kelvin

    // pulse parameters
    double pulse_start_time = 0.0; // Time when pulses begin in seconds
    double pulse_duration = 0.15; // Duration of each pulse in seconds (width)
    double pulse_period = 0.2; // Period of the pulses in seconds

    int n = num_cells + 1;   // number of nodes
    double dx = domain_length / num_cells;  // cell spacing

    // --- Allocate temperature arrays ---
    double* T     = new double[n];
    double* T_new = new double[n];

    // Initial condition: T(x,0) = 300 K
    double initial_temperature = 300.0; // Initial temperature in Kelvin at t=0
    for (int i = 0; i < n; ++i) {
        T[i] = initial_temperature;
    }

   for (int t_step = 0; t_step < num_time_steps; ++t_step) {
        double current_time = t_step * dt;

        // Determine if the thruster is on or off based on the pulse parameters
        bool thruster_on = false;
        if (current_time >= pulse_start_time) {
            double time_since_start = current_time - pulse_start_time;
            double time_in_cycle = fmod(time_since_start, pulse_period);
            if (time_in_cycle < pulse_duration) {
                thruster_on = true;
            }
        }

        // Update left boundary condition based on thruster state dT/dx at x=0
        if (thruster_on) {
            double left_boundary_on = (hg/k) * (T0g - T[0]); // Left boundary condition when thruster is on
            T_new[1] = T[0] + left_boundary_on * dx; // Apply Neumann BC at left boundary
        } else {
            double left_boundary_off = 0; // Left boundary condition when thruster is off
            T_new[1] = T[0] + left_boundary_off * dx; // Apply Neumann BC at left boundary
        }

        // Update right boundary condition dT/dx = 0 at x=L
        double right_boundary = 0; // Right boundary condition
        T_new[n] = T[n-1] + right_boundary * dx; // Apply Neumann BC at right boundary

        // Update interior nodes using FTCS scheme
        for (int i = 1; i < n; ++i) {
            T_new[i] = T[i] + dt * alpha *
                       (T[i-1] - 2.0*T[i] + T[i+1]) / (dx*dx);
        }

        for (int i = 0; i < n; ++i) {
            T[i] = T_new[i];
        }

}