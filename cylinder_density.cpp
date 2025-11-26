#include <iostream>

int main() {
    double z = 1; // cylinder height
    double r = 0.5; // cylinder radius

    double dz = 0.02;
    double dr = 0.02;

    int nz = 51; // number of nodes in z
    int nr = 26; // number of nodes in r

    // allocate memory for number density array
    double *nd = new double[nz*nr];
    double *nd_new = new double[nz*nr];
    double *I = new double[nz*nr];
    double *L = new double[nz*nr];

    // initialize number density array
    for (int i = 0; i < nz; i++) {        // loop over z
        for (int j = 0; j < nr; j++) {    // loop over r

            int idx = i*nr + j;           // flattened index

            nd[idx] = 0.0;                // initialize to zero
            nd_new[idx] = 0.0;            // initialize to zero
            L[idx] = 0.0;                  // initialize to zero

            double r_val = j * dr;          // physical radius at index j

            // inlet condition
            if (i==0 & r_val<=0.04) {   // corresponds to z = 0 and r = 0, 0.02, 0.04
                nd[idx] = 100;    // source
                nd_new[idx] = 100;
            }

            // outlet condition
            if (i == nz-1 & 0.1<=r_val & r_val<=0.2) { // corresponds to z = 1 and r = 0.1, 0.12, ..., 0.2
                nd[idx] = 0.0;      // sink
                nd_new[idx] = 0.0;
            }

            // initialize identity matrix
            if (i == j) {
                I[idx] = 1.0;
            } else {
                I[idx] = 0.0;
            }
        }
    }

    // crank-nicholson solver for unsteady diffusion equation in cylindrical coordinates
    double D = 0.1; // diffusion coefficient

    // Laplacian matrix
    for (int i = 0; i < nz; i++) {        // loop over z
        for (int j = 0; j < nr; j++) {    // loop over r

            int idx = i*nr + j;           // flattened index
            double f = nd[idx]; // center
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
            double d = 0.0;
            
            // top
            if (j != nr-1) {
                int t = j+1;
                int idx_t = i*nr + t;
                a = nd[idx_t]; // top
            }

            // left
            if (i != 0) {
                b = nd[idx-1]; // left
            }

            // bottom
            if (j != 0) {
                int l = j-1;
                int idx_l = i*nr + l;
                c = nd[idx_l]; // bottom
            }

            // right
            if (i != nz-1) {
                d = nd[idx+1]; // right
            }

            double r_val = j * dr;          // physical radius at index j
            // compute L matrix
            L[idx] = (a-c)/(2*dr*r_val) + (a-2*f+c)/(dr*dr) + (b-2*f+d)/(dz*dz);

        }
    }

    double dt = 0.01; // time step
    int nsteps = 5000; // number of time steps

    for (int k=0; k<nsteps; k++) { // time steps
        for (int i = 0; i < nz; i++) {        // loop over z
            for (int j = 0; j < nr; j++) {    // loop over r

                int idx = i*nr + j;           // flattened index

                // update number density using crank-nicholson
                nd_new[idx] = nd[idx] + D * dt * L[idx];
            }
        }

        // update nd for next time step
        for (int i = 0; i < nz; i++) {        // loop over z
            for (int j = 0; j < nr; j++) {    // loop over r

                int idx = i*nr + j;           // flattened index

                nd_new[idx] = nd[idx] + I[idx] + D * dt * L[idx];
            }
        }
    }

    return 0;
}
