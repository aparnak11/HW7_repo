#include <iostream>   // for screen output
#include <fstream>    // for file writing
#include <cmath>      // for sqrt, fmod

int main() {
    // simulation parameters
    double L = 0.01;   // Length of the domain in meters (x and y)
    double dt = 1e-4;              // Time step in seconds
    int num_cells = 100;           // Total spatial cells (we'll use sqrt for 2D)
    int num_time_steps = 15000;    // Number of time steps

    // material properties
    double K   = 65.0;     // Thermal conductivity [W/(m·K)]
    double alpha = 1.7e-5; // Thermal diffusivity [m^2/s]
    double hg  = 1.5e4;    // Convective heat transfer coefficient [W/(m^2·K)]
    double T0g = 2500.0;   // Gas stagnation temperature [K]

    // pulse parameters
    double pulse_start_time = 0.0;  // Time when pulses begin [s]
    double pulse_duration   = 0.15; // Duration of each pulse [s]
    double pulse_period     = 0.2;  // Period of the pulses [s]

    // 1D domain parameters
    int n = num_cells + 1;  // nodes in x direction
    double dx = L / num_cells;  // cell spacing in x

    // --- Allocate temperature matrix ---
    // Allocate an array of pointers (for rows)
    int rows = num_time_steps + 1;
    int cols = n;
    double** T = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        T[i] = new double[cols];
    }

    // Initial condition: T(t=0, x) = 300 K
    double initial_temperature = 300.0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            T[i][j] = initial_temperature;
        }
    }

    // ----- Main Loop -----
    for (int k = 0; k < num_time_steps; ++k) {
        // logic to determine if thruster is on or off
        double current_time = k * dt;

        bool thruster_on = false;

        if (current_time >= pulse_start_time) {
            double time_since_start = current_time - pulse_start_time;

            // Time within the current pulse cycle [0, pulse_period)
            double time_in_cycle = fmod(time_since_start, pulse_period);

            // Thruster is ON during the first 'pulse_duration' of each period
            thruster_on = (time_in_cycle < pulse_duration);
        }

        // spatial loop
        for (int i = 0; i < n; ++i) {
            // boundary conditions
            if (i == 0) { // left boundary
                if (thruster_on) {
                    double dTdx = (hg / K) * (T0g - T[k][i]);
                    T[k+1][0] = T[k][0] + dTdx * dx;
                    continue;
                } else {
                    T[k+1][i] = T[k][i]; // insulated boundary
                    continue;
                }
            }
            if (i == n-1) { // right boundary
                T[k+1][i] = T[k][i]; // insulated boundary
                continue;
            }

            // interior points
            T[k+1][i] = T[k][i] + alpha * dt * 
                (T[k][i+1] - 2*T[k][i] + T[k][i-1]) / (dx * dx);
        }

    }

    // ----- Output VTI file -----
    int nt_vis  = 100;                 // desired number of time samples for plotting
    int stride  = rows / nt_vis;
    if (stride < 1) stride = 1;
    nt_vis = rows / stride;            // recompute actual number of slices we’ll output

    double total_time = num_time_steps * dt;
    double dy_vis = (nt_vis > 1) ? total_time / (nt_vis - 1) : dt;  // spacing in time for the VTI


    std::ofstream out("field.vti");

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

    // x: space, 0..cols-1, spacing = dx
    // y: time (downsampled), 0..nt_vis-1, spacing = dy_vis
    // z: single slice
    out << "<ImageData WholeExtent=\"0 " << cols - 1
        << " 0 " << nt_vis - 1
        << " 0 0\" "
        << "Origin=\"0 0 0\" "
        << "Spacing=\"" << dx << " " << dy_vis << " 1\">\n";

    out << "<Piece Extent=\"0 " << cols - 1
        << " 0 " << nt_vis - 1
        << " 0 0\">\n";

    out << "<PointData Scalars=\"T\">\n";
    out << "<DataArray type=\"Float64\" Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    // ivis indexes the *visual* time slices
    for (int ivis = 0; ivis < nt_vis; ++ivis) {
        int k = ivis * stride;   // map to actual simulation time index (row in T)

        for (int i = 0; i < cols; ++i) {   // spatial index
            out << T[k][i] << " ";
        }
        out << "\n";
    }

    out << "</DataArray>\n";
    out << "</PointData>\n";
    out << "</Piece>\n";
    out << "</ImageData>\n";
    out << "</VTKFile>\n";

    out.close();

    // ----- Free memory -----
    for (int i = 0; i < rows; i++) {
        delete[] T[i];
    }
    delete[] T;
    return 0;
}
