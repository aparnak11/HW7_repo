/*
 * ASTE-499 Class 11 example
 * Unsteady advection-diffusion equation solver
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <chrono>
#include "Vec.h"
#include "Mat.h"

using namespace std;

struct World {
  double3 x0;	// origin
  double3 dh;   // cell spacing
  int3 nn;  	// number of nodes
  int U(int i, int j, int k=0) {return k*nn[1]*nn[0]+j*nn[0]+i;}
};

enum class AdScheme {UPWIND, CENTRAL};

//prototypes
void saveVTI(int ts, const World &world, const vector<double> &den,
			 const vector<double> &u, const vector<double> &v);
void solveGS(const Mat &mat, vector<double>&x, const vector<double>&b,
		     int max_it = 10000, double tol = 1e-4);


int main() {
 World world;
 world.x0 = {0,0,0};		// mesh origin
 world.dh = {0.2,0.2,0.1};	// cell spacing
 world.nn = {81,41,1};		// number of nodes in i, j, k
 double dt = 1e-2;

 int ni = world.nn[0];
 int nj = world.nn[1];
 int nk = world.nn[2];
 int nu = ni*nj*nk;

 double dx = world.dh[0];
 double dy = world.dh[1];

 SparseMat L(nu);		// diffusion operator

 vector<double> den(nu);
 vector<double> R(nu);

 vector<double> vel_u(nu); //velocity
 vector<double> vel_v(nu);

 const double D = 0.1;		// diffusion coefficient

 // set matrix coefficients for a 2D problem
 for (int i=0;i<ni;i++)
	 for (int j=0;j<nj;j++)
	 {
		 int u = world.U(i,j);
		 if (i == 0) {/* */}
		 else if (i==ni-1) {/* */}
		 else if (j==0) {/* */}
		 else if (j==nj-1) {/* */den[u]=0;}
		 else {
			 // standard 5-point stencil (2D)
			 double dx = world.dh[0];
			 double dy = world.dh[1];
			 L(u,u-ni) = D/(dy*dy);
			 L(u,u-1) = D/(dx*dx);
			 L(u,u) = (-2*D)*(1/(dx*dx) + 1/(dy*dy));
			 L(u,u+1) = D/(dx*dx);
			 L(u,u+ni) = D/(dy*dy);
		 }
	 }

  // add drop with density "10" at location i=0.25*ni, j=0.5*nj (quarter of the way in x, middle in y)
 den[world.U(0.25*ni,0.5*nj)] = 10.0;
 den[world.U(0.5*ni,0.1*nj)] = 10.0;


 // set flow velocity
 for (int i=0;i<ni;i++)
	 for (int j=0;j<nj;j++) {
		 int r = world.U(i,j);
		 //vel_u[r] = i/(ni-1.0)*5.0;
		 vel_u[r] = 1;
		 vel_v[r] = 0;
	 }

 auto t1 = chrono::high_resolution_clock::now();

  /*
  * (den_new-den)/dt + A*den = b
  * den_new - den = b-A*den*dt +den
  *
  */

 vector<double> div_vc(nu);

 // use dt=1e-3, D=0.1
 AdScheme ad_scheme = AdScheme::UPWIND;


 for (int ts=0;ts<1000;ts++) {
	 //compute div(v*c)
	 for (int i=1;i<ni-1;i++)
		 for (int j=1;j<nj-1;j++) {
			 int u = world.U(i,j);
			 double duc_dx, dvc_dy;

			 switch (ad_scheme) {
			 case AdScheme::UPWIND:
				 if (vel_u[u]>=0)
					 duc_dx = (vel_u[u]*den[u]-vel_u[u-1]*den[u-1])/dx;
				 else
					 duc_dx = (vel_u[u+1]*den[u+1]-vel_u[u]*den[u])/dx;
				 if (vel_v[u]>=0)
					 dvc_dy = (vel_v[u]*den[u]-vel_v[u-ni]*den[u-ni])/dy;
				 else
					 dvc_dy = (vel_v[u+ni]*den[u+ni]-vel_v[u]*den[u])/dy;
				 break;
			 case AdScheme::CENTRAL:
				 duc_dx = (vel_u[u+1]*den[u+1]-vel_u[u-1]*den[u-1])/(2*dx);
				 dvc_dy = (vel_v[u+ni]*den[u+ni]-vel_v[u-ni]*den[u-ni])/(2*dy);
				 break;
			 }



			 div_vc[u] = duc_dx+dvc_dy;
		 }


	 den = den + (R + L*den - div_vc)*dt;
	 if (ts%10==0) { cout<<ts<<endl;saveVTI(ts,world,den,vel_u,vel_v);}
 }
 // solves system A*den=b
 //solveGS(A,T,b, 10000, 1e-5);

 auto t2 = chrono::high_resolution_clock::now();

 chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
 //cout<<"Done in "<<time_span.count()<<" seconds using "<<A.memsize()/(1024.0*1024)<<" Mb of RAM"<<endl;
 return 0;
}

// solves A*x=b using GS-SOR method
void solveGS(const Mat &A, vector<double>&x, const vector<double>&b, int max_it, double tol) {
	cout<<"In solveGS"<<endl;
	for (int n=0;n<max_it;n++) {
		for (size_t r=0;r<A.nu;r++) {
			double dot = A.dot(x,r);  // A[r] times x dot product
			double g = (b[r] - (dot - A(r,r)*x[r]))/A(r,r);
			x[r] = x[r] + 1.4*(g-x[r]);		// SOR
		}

		if (n%25==0) {
			double sum = 0;
			for (size_t u=0;u<A.nu;u++) {
				double r = b[u]-A.dot(x,u);
				sum+=r*r;
			}
			double norm = sqrt(sum/A.nu);
			if (n%500==0) cout<<"n: "<<n<<", L2: "<<norm<<endl;
			if (norm<tol) {cout<<"Converged in "<<n<<" iterations"<<endl;break;}
		}
	}
}

// writes out density field to a "results" folder (needs to be manually created!)
void saveVTI(int ts, const World &world, const vector<double> &den,
		const vector<double> &u, const vector<double> &v)  {
	// output data
	stringstream ss;
	ss<<"results/field_"<<setfill('0')<<setw(6)<<ts<<".vti";

	ofstream out(ss.str());            //open output file
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData WholeExtent=\"0 "<<world.nn[0]-1<<" 0 "<<world.nn[1]-1<<" 0 "<<world.nn[2]-1<<"\" "
	   <<"Origin=\""<<world.x0<<"\" Spacing=\""<<world.dh<<"\">\n";
	out<<"<Piece Extent=\"0 "<<world.nn[0]-1<<" 0 "<<world.nn[1]-1<<" 0 "<<world.nn[2]-1<<"\">\n";

	out<<"<PointData>\n";
	
	out<<"<DataArray Name=\"den (#/m^3)\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">\n";
	for (const double &val:den) out<<val<<" ";
	out<<"</DataArray>\n";

	out<<"<DataArray Name=\"vel (#/m^3)\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
	for (size_t i=0;i<u.size();i++) out<<u[i]<<" "<<v[i]<<" 0 ";
	out<<"</DataArray>\n";

	out<<"</PointData>\n";

	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
}
