// Header file containing the definition of the structure of our Diff class, that we edit and further define in the diffraction.cpp file
// Linked by include-statement in program files that wishes to make use of this class to run simulations

#ifndef __diffraction_hpp__
#define __diffraction_hpp__
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class Diff{
    public:
    // we define a couple variables that'll help us discretize the space we're going to be working with 
    double h; // space step length for our xy-space
    double dt; // time step length
    double T; // total duration in time 
    int M; // variable defining the discretized space, equals the amount of points in our space
    // matrices needed to define our simulation space
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;
    arma::cx_mat U;
    arma::cx_cube S;     
    arma::mat V;    
    // vectors/matrices needed to further define our space
    arma::cx_vec u_1;
    arma::cx_vec u_0;
    arma::cx_double r; 
    
    Diff(double h, double dt, double T); // initializer, that will further define needed matrices and vectors, as well as a coupåle variables
    int indexes(int i, int j); // function defining indices, to navigate the space by i and j 
    void AB_def(); // function that defines the matrices A and B needed to perform the Crank-Nicolson Scheme
    void initialize(double xc, double sigmax, double px, double yc, double sigmay, double py); // 
    void pot_def(int slits, int V0); // function defining the potential that is used to model the wall the wave-function is going to interact with
    void simulate(); // function that simulates the evolution of our system with time
};