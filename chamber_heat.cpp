#include <iostream>
#include <fstream>          // for file writing

int main() {
    // simulation parameters
    double L = 0.01; // length in meters
    double dt = 1e-4; // time step in seconds
    double ncells = 100; // number of cells
    double nsteps = 15000; // number of time steps
    double dx = L / ncells; // cell size in meters

    // forward time with initial conditions
    double *T = new double[ncells*nsteps];
    double *Tnew = new double[ncells*nsteps];
    for (int i=0; i<ncells*nsteps; i++) {
        T[i] = 300.0; // initial temperature in K
        Tnew[i] = 300.0; // initial temperature in K
    }

    // pulse parameters
    double pulse_start_time = 0.0; // pulse start time in seconds
    double pulse_duration = 0.15; // pulse end time in seconds
    double pulse_period = 0.2; // pulse period in seconds

    // material properties
    double k_cond = 65.0; // thermal conductivity in W/m-K
    double alpha = 1.7e-5; // thermal diffusivity in m^2/s
    double hg = 1.5e4; // convective heat transfer coefficient in W/m^2-K
    double Tog = 2500.0; // hot gas temperature in K

    // boof pulse schedule
    int *x = new int[nsteps+1];
    for (int i=0; i<=nsteps; i++) {
        x[i] = 0; // initialize all time steps to thruster off
    }

    double *t = new double[nsteps+1];

    for (int k=0; k<=nsteps; k++) {
        t[k] = k * dt; // time array

        bool thruster_on = false;
        double current_time = k * dt;
        double time_in_cycle = fmod(current_time, pulse_period);
        if (time_in_cycle >= pulse_start_time && time_in_cycle < (pulse_start_time + pulse_duration)) {
            thruster_on = true; // thruster is on during this time step
        }

        // set pulse schedule
        if (thruster_on) {
            x[k] = 1;
        } else {
            x[k] = 0;
        }

        // main heat map loop
        for (int i=0; i<=ncells; i++) {
            if (i == 0) { // left boundary
                if (thruster_on) {
                    // thruster on boundary condition
                    double dTdx = (hg/k_cond) * (Tog - T[i]);
                    Tnew[i] = T[i] +  alpha * dt * (T[i+1] - 2*T[i]) * dTdx / dx; // forward difference
                    continue;
                } else {
                    // thruster off boundary condition
                    Tnew[i] = T[i];
                    continue;
                }
            } else if (i == ncells) { // right boundary
                Tnew[i] = T[i]; // insulated boundary condition
                continue;
            }

            // interior points
            Tnew[i] = T[i] + alpha * dt * (T[i+1] - 2*T[i] + T[i-1]) / (dx * dx);

            //T[i] = Tnew[i]; // update temperature array for next time step
        }
    }

    // vti output
    double x0 = 0.1;
    double y0 = 0.1;

    std::ofstream out("field.vti");

    out<<"<VTKFile type=\"ImageData\">\n";
    out << "<ImageData WholeExtent=\"0 " << ncells << " 0 " << ncells << " 0 0\"";
    out<<" Origin=\""<<x0<<" "<<y0<<" "<<0.0<<"\"";
    out<<" Spacing=\""<<dx<<" " <<dt<<" "<<0.0<<"\">\n";
    out<<"<Piece Extent=\"0 "<<ncells<<" 0 "<<nsteps<<" 0 "<<0<<"\">\n";
    out<<"<PointData>\n";

    out<<"<DataArray Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
    for (int n=0; n<ncells*nsteps; n++) out<<T[n]<<" ";
    out<<"\n</DataArray>\n";

    out<<"</PointData>\n";
    out<<"</Piece>\n";
    out<<"</ImageData>\n";
    out<<"</VTKFile>\n";

    // ---- Dump pulse schedule to CSV ---- //
    std::ofstream file("pulse.csv");
    file << "time,thruster\n";   // header

    for (int k=0; k<=nsteps; k++) {
        file << t[k] << "," << x[k] << "\n";
    }
    file.close();

    // free memory
    delete[] t;
    delete[] x;
    delete[] T;
    delete[] Tnew;


    return 0;
}