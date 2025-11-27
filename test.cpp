#include <iostream>
#include <fstream> // for file writing

// simulation parameters
const double L = 0.01;        // length in meters
const double dt = 1e-4;       // time step in seconds
const int ncells = 101;       // number of cells
const int nsteps = 15001;     // number of time steps
const double dx = L / ncells; // cell size in meters

// pulse parameters
const double pulse_start_time = 0.0; // pulse start time in seconds
const double pulse_duration = 0.15;  // pulse end time in seconds
const double pulse_period = 0.2;     // pulse period in seconds

// material properties
const double k_cond = 65.0;  // thermal conductivity in W/m-K
const double alpha = 1.7e-5; // thermal diffusivity in m^2/s
const double hg = 1.5e4;     // convective heat transfer coefficient in W/m^2-K
const double Tog = 2500.0;   // hot gas temperature in K

double A[ncells][nsteps];
double T[ncells * nsteps];

bool is_thruster_on(int k, double dt)
{

    bool thruster_on = false;
    double current_time = k * dt;
    double time_in_cycle = fmod(current_time, pulse_period);
    if (time_in_cycle >= pulse_start_time && time_in_cycle < (pulse_start_time + pulse_duration))
    {
        thruster_on = true; // thruster is on during this time step
    }
    return thruster_on;
}

int main()
{
    std::cout << "Starting chamber heat simulation..." << std::endl;

    for (int i = 0; i < ncells; i++)
    {
        for (int k = 0; k < nsteps; k++)
        {
            A[i][k] = 300.0; // initial temperature in K
        }
    }
    std::cout << "Initialized temperature array." << std::endl;

    for (int k = 1; k < nsteps; k++)
    { // time steps
        for (int i = 0; i < ncells; i++)
        { // spatial loop
            // boundary conditions
            // left boundary
            if (i == 0)
            {
                if (is_thruster_on(k, dt))
                {
                    double dTdx = (hg / k_cond) * (Tog - A[i][k - 1]);
                    double l_ghost = A[i + 1][k - 1] + 2 * dx * dTdx; // ghost cell
                    A[i][k] = A[i][k - 1] + alpha * dt * (l_ghost - 2 * A[i][k - 1] + A[i + 1][k - 1]) / (dx * dx);
                    if (k < 2000)
                    {
                        std::cout << "Thruster ON at time step " << k << ", cell 0" << A[i][k] << std::endl;
                    }
                }
                else
                {
                    // thruster off boundary condition
                    double l_ghost = A[i + 1][k - 1]; // insulated boundary condition
                    A[i][k] = A[i][k - 1] + alpha * dt * (l_ghost - 2 * A[i][k - 1] + A[i + 1][k - 1]) / (dx * dx);
                    if (k < 2000)
                    {
                        std::cout << "Thruster OFF at time step " << k << ", cell 0" << A[i][k] << std::endl;
                    }
                }
            }
            // right boundary
            else if (i == ncells - 1)
            {
                // double r_ghost = A[i - 1][k - 1]; // insulated boundary condition
                double r_ghost = 300.0;
                A[i][k] = A[i][k - 1] + alpha * dt * (A[i - 1][k - 1] - 2 * A[i][k - 1] + r_ghost) / (dx * dx);
            }
            else
            {
                // interior cells
                A[i][k] = A[i][k - 1] + alpha * dt * (A[i - 1][k - 1] - 2 * A[i][k - 1] + A[i + 1][k - 1]) / (dx * dx);
            }
        }
    }
    std::cout << "Completed time stepping." << std::endl;

    std::cout << "A[0][0] = " << A[0][0]
              << ", A[0][1000] = " << A[0][1000]
              << ", A[100][15000] = " << A[100][15000]
              << ", A[50][1000] = " << A[50][1000] << std::endl;

    for (int k = 0; k < nsteps; k++)
    { // time steps
        for (int i = 0; i < ncells; i++)
        { // spatial loop
            int idx = k * ncells + i;
            T[idx] = A[i][k];
        }
    }
    std::cout << "Prepared data for output." << std::endl;

    // ==================== DOWNSAMPLE TO EVERY 100th TIMESTEP ====================

    // Keep every 100th timestep
    const int stride = 100;
    const int nsteps_out = nsteps / stride; // → 150
    std::vector<double> T_out(ncells * nsteps_out);

    // Flatten A[][] into T_out[] only for selected timesteps
    int idx = 0;
    for (int k = 0; k < nsteps; k += stride) // sample time dimension
    {
        for (int i = 0; i < ncells; i++) // full spatial dimension
        {
            T_out[idx++] = A[i][k];
        }
    }

    std::cout << "Downsampled to " << nsteps_out << " timesteps." << std::endl;

    // ============================== WRITE VTI ===================================

    double x0 = 0.0, y0 = 0.0;
    std::ofstream out("field_downsampled.vti");

    out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

    out << "<ImageData WholeExtent=\"0 " << ncells - 1
        << " 0 " << nsteps_out - 1
        << " 0 0\" Origin=\"" << x0 << " " << y0 << " 0\" "
        << "Spacing=\"" << dx << " " << dt * stride << " 0\">\n";

    out << "<Piece Extent=\"0 " << ncells - 1
        << " 0 " << nsteps_out - 1
        << " 0 0\">\n";

    out << "<PointData Scalars=\"T\">\n<DataArray Name=\"Temperature\" "
        << "NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";

    for (double v : T_out)
        out << v << " ";

    out << "\n</DataArray>\n</PointData>\n<CellData/>\n</Piece>\n</ImageData>\n</VTKFile>\n";

    std::cout << "VTK export complete → field_downsampled.vti" << std::endl;

    return 0;
}
