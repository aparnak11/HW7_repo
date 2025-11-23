#include <iostream>   // for screen output
#include <fstream>    // for file writing
#include <cmath>      // for sqrt, fmod

int main() {
    double domain_length = 0.01;   // Length of the domain in meters (x and y)
    double dt = 1e-4;              // Time step in seconds
    int num_cells = 100;           // Total spatial cells (we'll use sqrt for 2D)
    int num_time_steps = 15000;    // Number of time steps

    // material properties
    double k   = 65.0;     // Thermal conductivity [W/(m·K)]
    double alpha = 1.7e-5; // Thermal diffusivity [m^2/s]
    double hg  = 1.5e4;    // Convective heat transfer coefficient [W/(m^2·K)]
    double T0g = 2500.0;   // Gas stagnation temperature [K]

    // pulse parameters
    double pulse_start_time = 0.0;  // Time when pulses begin [s]
    double pulse_duration   = 0.15; // Duration of each pulse [s]
    double pulse_period     = 0.2;  // Period of the pulses [s]

    // 2D grid: assume num_cells = (ni-1)*(nj-1) and square grid
    int ni = static_cast<int>(std::sqrt(num_cells)) + 1;  // nodes in x
    int nj = static_cast<int>(std::sqrt(num_cells)) + 1;  // nodes in y

    double x0 = 0.1;  // origin x (for VTK)
    double y0 = 0.1;  // origin y

    double dx = domain_length / (ni - 1);  // cell spacing in x
    double dy = domain_length / (nj - 1);  // cell spacing in y

    int nn = ni * nj;  // total number of nodes

    // --- Allocate temperature arrays ---
    double* T     = new double[nn];
    double* T_new = new double[nn];

    // Initial condition: T(x,y,0) = 300 K
    double initial_temperature = 300.0;
    for (int n = 0; n < nn; ++n) {
        T[n] = initial_temperature;
    }

    double dxsqr = dx * dx;
    double dysqr = dy * dy;

    // ----- Time marching -----
    for (int t_step = 0; t_step < num_time_steps; ++t_step) {
        double current_time = t_step * dt;

        // Determine if the thruster is on or off based on the pulse parameters
        bool thruster_on = false;
        if (current_time >= pulse_start_time) {
            double time_since_start = current_time - pulse_start_time;
            double time_in_cycle = std::fmod(time_since_start, pulse_period);
            if (time_in_cycle < pulse_duration) {
                thruster_on = true;
            }
        }

        // --------------------------------------------------------
        // 1) Apply boundary conditions on T (OLD field)
        // --------------------------------------------------------

        // Left boundary (i = 0): convective BC when thruster ON, zero-gradient when OFF
        for (int j = 0; j < nj; ++j) {
            int i0 = 0;
            int i1 = 1;
            int id0 = j * ni + i0;
            int id1 = j * ni + i1;

            if (thruster_on) {
                // dT/dx = (hg/k)*(T0g - T(i=0))
                double dTdx = (hg / k) * (T0g - T[id0]);
                // (T1 - T0)/dx = dTdx  ->  T1 = T0 + dx * dTdx
                T[id1] = T[id0] + dx * dTdx;
            } else {
                // dT/dx = 0  ->  T1 = T0
                T[id1] = T[id0];
            }
        }

        // Right boundary (i = ni-1): zero-gradient dT/dx = 0
        for (int j = 0; j < nj; ++j) {
            int iR  = ni - 1;
            int iRm = ni - 2;
            int idR  = j * ni + iR;
            int idRm = j * ni + iRm;
            T[idR] = T[idRm];
        }

        // Bottom boundary (j = 0): zero-gradient dT/dy = 0
        for (int i = 0; i < ni; ++i) {
            int j0 = 0;
            int j1 = 1;
            int id0 = j0 * ni + i;
            int id1 = j1 * ni + i;
            T[id0] = T[id1];
        }

        // Top boundary (j = nj-1): zero-gradient dT/dy = 0
        for (int i = 0; i < ni; ++i) {
            int jT  = nj - 1;
            int jTm = nj - 2;
            int idT  = jT * ni + i;
            int idTm = jTm * ni + i;
            T[idT] = T[idTm];
        }

        // --------------------------------------------------------
        // 2) FTCS update for interior nodes (1 <= i <= ni-2, 1 <= j <= nj-2)
        // --------------------------------------------------------
        for (int j = 1; j < nj - 1; ++j) {
            for (int i = 1; i < ni - 1; ++i) {
                int id  = j * ni + i;
                int idL = j * ni + (i - 1);
                int idR = j * ni + (i + 1);
                int idD = (j - 1) * ni + i;
                int idU = (j + 1) * ni + i;

                double d2Tdx2 = (T[idL] - 2.0 * T[id] + T[idR]) / dxsqr;
                double d2Tdy2 = (T[idD] - 2.0 * T[id] + T[idU]) / dysqr;

                T_new[id] = T[id] + dt * alpha * (d2Tdx2 + d2Tdy2);
            }
        }

        // --------------------------------------------------------
        // 3) Copy boundary values into T_new directly from T
        // (they already satisfy the BCs)
        // --------------------------------------------------------
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                if (i == 0 || i == ni - 1 || j == 0 || j == nj - 1) {
                    int id = j * ni + i;
                    T_new[id] = T[id];
                }
            }
        }

        // --------------------------------------------------------
        // 4) Swap T_new → T for the next time step
        // --------------------------------------------------------
        for (int n = 0; n < nn; ++n) {
            T[n] = T_new[n];
        }
    }

    // ----- Output VTI file -----
    std::ofstream out("field.vti");

    out << "<VTKFile type=\"ImageData\">\n";
    out << "<ImageData WholeExtent=\"0 " << ni-1
        << " 0 " << nj-1 << " 0 0\"";
    out << " Origin=\"" << x0 << " " << y0 << " " << 0.0 << "\"";
    out << " Spacing=\"" << dx << " " << dy << " " << 0.0 << "\">\n";
    out << "<Piece Extent=\"0 " << ni-1
        << " 0 " << nj-1 << " 0 0\">\n";
    out << "<PointData>\n";

    out << "<DataArray Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
    for (int n = 0; n < nn; n++) out << T[n] << " ";
    out << "\n</DataArray>\n";

    out << "</PointData>\n";
    out << "</Piece>\n";
    out << "</ImageData>\n";
    out << "</VTKFile>\n";

    // free arrays
    delete[] T;
    delete[] T_new;

    return 0; // normal exit
}
