#include <iostream>
#include <fstream>
#include <cmath>

#define IDX(n,i,N) ((n)*(N) + (i))   // flatten time n, space i

int main() {

    // Parameters
    double L  = 0.01;
    int    N  = 101;
    double dx = L / (N - 1);
    double dt = 1e-4;
    int    nt = 15000;

    double alpha = 1.7e-5;
    double k     = 65.0;
    double hg    = 1.5e4;
    double T0g   = 2500.0;

    double r = alpha * dt / (dx * dx);

    // allocate memory: full 2D field (time × space)
    double *T    = new double[N * nt];
    double *Tnew = new double[N];

    // initial condition at t = 0
    for (int i = 0; i < N; i++)
        T[IDX(0,i,N)] = 300.0;

    // ---- TIME LOOP ----
    for (int n = 0; n < nt-1; n++) {

        double time = n * dt;
        bool thruster_on = fmod(time, 0.2) < 0.15;

        double* Told  = &T[IDX(n,0,N)];
        double* Tnext = Tnew;

        // ghost node at x<0
        double Tghost;
        if (thruster_on) {
            Tghost = Told[1] - 2*dx*(hg/k)*(T0g - Told[0]);
        } else {
            Tghost = Told[1];
        }

        // left boundary
        Tnext[0] = Told[0] + r*(Told[1] - 2*Told[0] + Tghost);

        // interior
        for (int i = 1; i < N-1; i++) {
            Tnext[i] = Told[i] + r*(Told[i+1] - 2*Told[i] + Told[i-1]);
        }

        // right boundary (Neumann 0)
        Tnext[N-1] = Tnext[N-2];

        // store next row
        for (int i = 0; i < N; i++)
            T[IDX(n+1,i,N)] = Tnext[i];
    }

    // total simulated time (for reference)
    double total_time = (nt - 1) * dt;   // ~1.5 s

    // ---- WRITE VTI ----
    std::ofstream out("field.vti");

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

    // WholeExtent: index ranges
    // x: 0..N-1 (space), y: 0..nt-1 (time), z: 0..0
    out << "<ImageData WholeExtent=\"0 " << (N-1)
        << " 0 " << (nt-1)
        << " 0 0\" "
        // Origin: (x=0 m, t=0 s, z=0)
        << "Origin=\"0 0 0\" "
        // Spacing: dx in meters, dt in seconds, 1 in z
        << "Spacing=\"" << dx << " " << dt << " 1\">\n";

    out << "<Piece Extent=\"0 " << (N-1)
        << " 0 " << (nt-1)
        << " 0 0\">\n";

    out << "<PointData Scalars=\"Temperature\">\n";
    out << "<DataArray type=\"Float64\" Name=\"T\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    // row-major: loop over time (rows), then space (columns)
    for (int n = 0; n < nt; n++) {
        for (int i = 0; i < N; i++)
            out << T[IDX(n,i,N)] << " ";
        out << "\n";
    }

    out << "</DataArray>\n";
    out << "</PointData>\n";
    out << "</Piece>\n";
    out << "</ImageData>\n";
    out << "</VTKFile>\n";

    out.close();

    delete[] T;
    delete[] Tnew;

    return 0;
}