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

    double initial_temperature = 300.0; // Initial temperature in Kelvin at t=0
    // boundary conditions
    double left_boundary_on = (hg/k) * (T0g - T_t0); // Left boundary condition when thruster is on
    double left_boundary_off = 0; // Left boundary condition when thruster is off
    double right_boundary = 0; // Right boundary condition

    // pulse parameters
    double pulse_start_time = 0.0; // Time when pulses begin in seconds
    double pulse_duration = 0.15; // Duration of each pulse in seconds (width)
    double pulse_period = 0.2; // Period of the pulses in seconds

}