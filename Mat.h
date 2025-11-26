#ifndef _MAT_H
#define _MAT_H

#include <vector>
#include <stdexcept>
#include <cassert>

#include "Vec.h"

// vector subtraction c = a-b
std::vector<double> operator-(const std::vector<double>&a, const std::vector<double>&b) {
	std::vector<double> c(a.size());
	for (size_t i=0;i<a.size();i++) c[i] = a[i]-b[i];
	return c;
}

// vector addition c = a+b
std::vector<double> operator+(const std::vector<double>&a, const std::vector<double>&b) {
	std::vector<double> c(a.size());
	for (size_t i=0;i<a.size();i++) c[i] = a[i]+b[i];
	return c;
}

// vector scaling c = a*s
std::vector<double> operator*(const std::vector<double>&a, double s) {
	std::vector<double> c(a.size());
	for (size_t i=0;i<a.size();i++) c[i] = a[i]*s;
	return c;
}

// matrix base class, specifies the interface to be implemented
class Mat {
public: 
  Mat(size_t nu) : nu{nu} {}
  virtual ~Mat() {}
  virtual double& operator()(int i, int j) = 0;
  virtual double operator()(int i, int j) const = 0;
  virtual double dot(const std::vector<double> &x, size_t r) const = 0;
  virtual size_t memsize() = 0;  			// matrix size in bytes
  virtual void clearRow(int r) = 0;     // sets all coeffs of row r to 0
  const size_t nu;		// number of rows (unknowns)

  //converts the matrix to an identity matrix
  void makeIdentity() {
    for (int r=0;r<nu;r++) {clearRow(r); operator()(r,r) = 1.0;}
  }

// computes b=A*x
  friend std::vector<double> operator* (const Mat &A, const std::vector<double> &x) {
    std::vector<double> b(A.nu);
    for (size_t u=0;u<A.nu;u++)
      b[u] = A.dot(x,u);
    return b;
  }

};


/* dense matrix, only good for small matrixes!*/
class DenseMat : public Mat {
public:
  DenseMat(size_t nu) : Mat{nu} {
    // memory for coefficients
    a = new double*[nu];
    for (size_t i=0;i<nu;i++)
    a[i] = new double[nu];

    // initialize all to zero
    for (size_t i=0;i<nu;i++)
      for (size_t j=0;j<nu;j++)
        a[i][j] = 0;
    }

  ~DenseMat() {
    for (size_t i=0;i<nu;i++)
      delete[] a[i];
    delete[] a;
  }

  // returns reference to the A[i,j] coefficient in the full matrix
  double& operator()(int i, int j) {return a[i][j];}
  double operator()(int i, int j) const {return a[i][j];}

  // returns dot product of matrix row r with vector x
  double dot(const std::vector<double> &x, size_t r) const {
    double sum = 0;
    for (size_t c=0;c<nu;c++)
      sum+=a[r][c]*x[c];
    return sum;
    }
	
    void clearRow(int r) {
      for (unsigned c=0;c<nu;c++) a[r][c] = 0;
    }
  
   // return memory size in bytes
   size_t memsize() {return nu*nu*sizeof(double);}

protected:
  double **a;		//coefficients, [nu][nu]
};


/* sparse matrix with up to 7 non-zero values per row*/
class SparseMat : public Mat {
public:

	SparseMat(size_t nu) : Mat{nu} {

		// memory for coefficients
		a = new double*[nu];
		for (size_t i=0;i<nu;i++)
			a[i] = new double[max_vals];

        // memory for column indexes
		c = new int*[nu];
			for (size_t i=0;i<nu;i++)
			c[i] = new int[max_vals];

		// clear data, column set to -1 are not set
		for (size_t i=0;i<nu;i++)
			for (size_t j=0;j<max_vals;j++) {
				c[i][j] = -1;  
			    a[i][j] = 0;   
			}
	}

	~SparseMat() {
		for (size_t i=0;i<nu;i++)
			delete[] a[i];
		delete[] a;

		for (size_t i=0;i<nu;i++)
			delete[] c[i];
		delete[] c;
		
	}

	// returns reference to the A[i,j] coefficient in the full matrix
	double& operator()(int i, int j) {
		//search for the sparse column corresponding to full matrix column j
		for (int v=0;v<max_vals;v++) {
			//did we reach an empty slot? If so, make it correspond to column j
			if (c[i][v]<0) c[i][v] = j;

			// does this sparse column map to j? If so, return the coeff
			if (c[i][v]==j) return a[i][v];
		}

		//getting here implies that all max_val slots are already occupied by columns other than j
		std::runtime_error("Sparse matrix too small!");
		return a[i][0];  // just to avoid a warning
	}

	// identical to the function above but for read-only access
	double operator()(int i, int j) const {
		//search for the sparse column corresponding to full matrix column j
		for (size_t v=0;v<max_vals;v++) {
			//did we reach an empty slot? If so, make it correspond to column j
			if (c[i][v]<0) c[i][v] = j;

			// does this sparse column map to j? If so, return the coeff
			if (c[i][v]==j) return a[i][v];
		}

		//getting here implies that all max_val slots are already occupied by columns other than j
		std::runtime_error("Sparse matrix too small!");
		return 0;
	}
  
    void clearRow(int r) {
	  for (int v=0;v<max_vals;v++) {a[r][v] = 0; c[r][v] = -1;}
	}

	// returns dot product of matrix row r with vector x
	double dot(const std::vector<double> &x, size_t r) const {
		double sum = 0;
		// loop up to max_vals time and until c[r][v] becomes negative
		for (size_t v=0;v<max_vals && c[r][v]>=0;v++)
			sum+=a[r][v]*x[c[r][v]];		// c[v] is effective the "j" in full matrix
		return sum;
	}

    friend SparseMat operator*(double s, const SparseMat &M) {
      size_t nu = M.nu;
      SparseMat A(nu);
      for (int r=0;r<nu;r++) 
        for (int c=0;c<nu;c++) A(r,c) = s*M(r,c);
      return A;
    }
    
    // matrix subtraction
    friend SparseMat operator-(const SparseMat &A, const SparseMat &B) {
      size_t nu = A.nu;
      SparseMat M(nu);
      
      // first copy A to M
      M = A;
      
      // now subtract B
      for (int r=0;r<nu;r++) {
        for (int v=0;v<max_vals;v++) {
          int c = B.c[r][v];
          if (c>=0) M(r,c) -= B.a[r][v];
          }
        }
       return M;
    }
    
    friend SparseMat operator+(const SparseMat &A, const SparseMat &B) {
      size_t nu = A.nu;
      SparseMat M(nu);
      
      // first copy A to M
      M = A;
      
      // now subtract B
      for (int r=0;r<nu;r++) {
        for (int v=0;v<max_vals;v++) {
          int c = B.c[r][v];
          if (c>=0) M(r,c) += B.a[r][v];
          }
        }
       return M;
    }
    
    void operator=(const SparseMat &A) {
      assert(nu==A.nu);
      for (int r=0;r<nu;r++)
        for (int v=0;v<max_vals;v++) {
          this->a[r][v] = A.a[r][v];
          this->c[r][v] = A.c[r][v];
        }
    }
    
   // return memory size in bytes
   size_t memsize() {return nu*max_vals*sizeof(double);}

protected:
    static constexpr int max_vals = 7;   //max non-zeros per row
    double **a;		//coefficients, [nu][7]
    int **c;        //columns in full matrix, [nu][7], -1 if not set
};



#endif
