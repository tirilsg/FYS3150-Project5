// Source file containing the further definition of the structure of our Diff class from the header file diffraction.hpp

#include "diffraction.hpp" // includes the definition we made in the header file 
#include <armadillo>
#include <iostream>
#include <complex>
#include <fstream>
#include <assert.h>

Diff::Diff(double h, double dt, double T) : h(h), dt(dt), T(T) {
    // defines the variables h, dt and T, as well as the amounts of steps in the xyspace M
    M = 1. / h + 1;
    r = arma::cx_double(0, dt / (2 * h * h));
    // defines the size of arrays that are to be containing the potential V, wave function U and the simulation data S
    U.zeros(M - 2, M - 2);
    V.zeros(M - 2, M - 2);
    S.zeros(M - 2, M - 2, T / dt + 1);
    // A and B, defined as cubes in the (M-2)*(M-2) x (M-2)*(M-2) space
    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

int Diff::indexes(int i, int j){ // definition of function translating the indices (i,j) to a single index k
  return i + (M - 2) * j;
}

void Diff::pot_def(int slits, int V0){// definition of potential as a function of amounts of slits
  double length = 0.05; //length of wall seperating slits
  double aperture = 0.05; // slit opening size in z-direction
  double thick = 0.02; // thickness of the wall in the x-dimension
  // definition of variables needed to further define areas where we wish to place a potential of a large size, so that the wave-function gets reflected where needed
  int cy = (M - 2) / 2; // centre point in the y-direction
  int cl = static_cast<int>(length / h); // length per step in space
  int a = static_cast<int>(aperture / h); // aperture per step in space
  int t = static_cast<int>(thick / h);   // thickness per step in space
  int nc = slits - 1; // implementing a variable for slits adjusted to our loop-definition
  int el = (M - 2) - slits * a - nc * cl;
  el /= 2;
  for (int i = 0; i < (M - 2); i++){ // i dimension loop
    for (int j = cy - t / 2; j <= cy + t / 2; j++){ // j dimension loop
      if (slits == 1){ // checks if the wall is to contain a singular slit
        if (i < el || i > el + a){ // if conditions are met
          V(i, j) = V0; // the potential at the coordinates (i,j) is set to V0
        }
      }
      else if (slits == 2){// checks if the wall is to contain double slits
        if (i < el || (i >= el + a && i <= el + a + cl) || i > el + 2 * a + cl){ // if conditions are met
          V(i, j) = V0; // the potential at the coordinates (i,j) is set to V0
        }
      }
      else if (slits == 3){// checks if the wall is to contain tripple slits
        if (i < el || (i >= el + a && i <= el + a + cl) ||(i > el + 2 * a + cl && i <= el + 2 * (a + cl)) || i > el + 2 * a + 3 * cl){
          V(i, j) = V0;
        }
      }
    }
  }
}

void Diff::AB_def(){ // definition of A and B vectors as defined in the report
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros); // definition of the shapes of vectors a and b
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros); // definition of a diagonal matrice 
  sub_diag.diag(-1).fill(r); // fills the matrice sub_diag with r in the subdiagonals
  sub_diag.diag(1).fill(r);
  A.diag(M - 2).fill(-r); // defines the diagonals of the matrices with -r for A
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r); // defines the diagonals of the matrices with r for B
  B.diag(2 - M).fill(r);
  arma::cx_double a_k, b_k; //defines the points a_k and b_k
  int k; // defines k 
  // loop defining the matrices a and b
  for (int i = 0; i < M - 2; i++){ // loops through points in the i dimension 
    for (int j = 0; j < M - 2; j++){ // loops through the points in j dimension 
      k = indexes(i, j); //defines k for this 
      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j)); // definition of  a_k as described in the report
      b_k = std::conj(a_k); // defines b_k as defined in the report
      a(k) = a_k; // sets a_k and b_k as the point in a and b matrices
      b(k) = b_k;
    }
  }
  for (int i = 0; i < M - 2; i++){ // another loop that defines the positions within the A and B matrices as the diagonal matrices
    int j = i * (M - 2);
    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag; // -diagonal matrice for A matrice
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag; // diagonal matrice for B matrice
  }
  A.diag() = a; // defines the diagonal of the A matrice by the matrice a
  B.diag() = b; // defines the diagonal of the B matrice by the matrice b
}

void Diff::simulate(){ // function that runs our simulations 
  for (int i = 1; i < (T / dt + 1); i++){
    u_0 = U.as_col(); // the old wave function that we wish to evolve a step
    u_1 = arma::spsolve(A, B * u_0); // new wave function that has been evolved
    for (int n = 0; n < U.n_cols; n ++){ // looping through columns
      for (int k = 0; k < M - 2; k++){ // looping through the matrice
        U.col(n)(k) = u_1(k + n * (M - 2)); // defines U by the evolved wave function estimation
      }
    }
    S.slice(i) = U; // stores the estimation in the cube S, that can be extracted later
    u_0 = u_1; // replaces the old definition with the new
  }
}

void Diff::initialize(double xc, double sigmax, double px, double yc, double sigmay, double py){ // function initializing our system
  double re, im; // defines doubles real and imaginary variables
  int xi, yi;
  for (int i = 0; i < M - 2; i++){ // i dimension loop 
    yi = i;
    for (int j = 0; j < M - 2; j++){ // j dimension loop
      xi = j;
      re = -(xi * h - xc) * (xi * h * xc) / (2 * sigmax * sigmax)-(yi * h - yc) * (yi * h - yc) / (2 * sigmay * sigmay); // Gaussian wave packet definition
      im = px * (xi * h - xc) + (yi * h - yc); // Gaussian wave packet definition
      U(i, j) = std::exp(arma::cx_double(re, im)); // saves the wave function in the indices (i,j)
    }
  }
  U /= std::sqrt(arma::accu(arma::conj(U) % U)); // normalizing 
  S.slice(0) = U; // sets initial condition at space 0 in ourcube
}
