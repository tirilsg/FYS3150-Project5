#ifndef __diffraction_hpp__
#define __diffraction_hpp__
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class Diff{
public:
    double h, dt, T;
    int M, N, slits;
    int V0; 
    arma::cx_double r;       
    arma::sp_cx_mat A, B;    
    arma::cx_mat U;         
    arma::cx_vec u_new, u_current;
    arma::mat V;           
    arma::cx_cube S;          
    Diff(double h_in, double dt_in, double T_in);
    int index_k(int i, int j);
    void fill_matrices();
    void simulate();
    void set_initial_state(double xc, double sigmax, double px, double yc, double sigmay, double py);
    void set_potential(int slits, int V0);
};











/*class Diff{
public:
    double h, dt, T;
    int M, N, slits;
    int V0; 
    arma::cx_double r;       
    arma::sp_cx_mat A, B;    
    arma::cx_mat U;         
    arma::cx_vec u_new, u_current;
    arma::mat V;           
    arma::cx_cube S;          
    Diff(double h_in, double dt_in, double T_in);
    int index_k(int i, int j);
    void fill_matrices();
    void simulate();
    void set_initial_state(double xc, double sigmax, double px, double yc, double sigmay, double py);
    void set_potential(int slits, int V0);
};*/




//////////////////////////////////////////////////////NOTES
/*class Diff{
public:
    double h, dt, T;
    int M, N, n_slits; 
    arma::cx_double r;       
    arma::sp_cx_mat A, B;    
    arma::cx_mat U;         
    arma::cx_vec u_new, u_current;
    arma::mat V;           
    arma::cx_cube S;          
    Diff(double h_in, double dt_in, double T_in);
    int index_k(int i, int j);
    void fill_matrices();
    void solve();
    void set_initial_state(double x_c, double sigma_x, double p_x, double y_c, double sigma_y, double p_y);
    void set_potential(int n_slits);
};*/


/*class Diff {
private:
    double r, dt, L, W, V0, b, h, m, x0, sigma, x_max, y_max, t_max, sigma_t;
    double n_slit, thickx, centerx, slit_sep, aperture, xc, yc, widthx, widthy, px, py;
    int M;
    arma::mat V;
    arma::mat u;
    arma::cx_double r;
    arma::mat A;
    arma::mat B;
    int len;
    double T;
    int n_timesteps;
    arma::cube storage;
    arma::cube Re;
    arma::cube Im;

public:
    Diff(double r_, double dt_, double L_, double W_, double V0_, double b_, double h_, double m_, double x0_, double sigma_, int M_, double x_max_, double y_max_, double t_max_, double sigma_t_,
         double n_slit_, double thickx_, double centerx_, double slit_sep_, double aperture_, double xc_, double yc_, double widthx_, double widthy_, double px_, double py_);

    arma::mat potential(double potential, int M, double n_slit, double thickx, double centerx, double slit_sep, double aperture);

    arma::mat wave_init(double xc, double yc, double widthx, double widthy, double px, double py, int M);

    void create_AB(const arma::mat& V, double h, double dt, int M);

    // Add other member functions as needed
};*/
///////////////////herfrsa
/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;
        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;
        arma::cube Re;
        arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        //void create_AB(arma::mat& V, const double h, const double dt, const int M);
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 
        double h;               
        int N;                  
        arma::cx_double r;         
        arma::sp_cx_mat A, B;          
        arma::cx_mat U;             
        arma::cx_vec u_new, u_current; 
        arma::cx_cube S;

        int index_k(int i, int j);                
};*/
/*class Matrix
{
public:
    // Member variables
    double h, dt, T;                // Spatial and time step length and total time
    int M, N;                       // No. of points along x, y and time axes

    arma::cx_double r;              // Predefined constant, r = i∆t/h^2

    arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
    arma::cx_mat U;                 // State matrix
    arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
    arma::mat V;                    // Potential matrix

    arma::cx_cube S;                // Storing states

    // Constructor
    Matrix(double h_in, double dt_in, double T_in, double x_c, double sigma_x, double p_x, double y_c,
           double sigma_y, double p_y, int n_slits);

    // Takes indices (i, j) and returns the corresponding k index
    int index_k(int i, int j);

    // Fill matrices A and B
    void fill_matrices();

    // Solve matrix equation Au^(n+1) = Bu^n
    void solve();

    // Set the initial state of the system
    void set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                           double sigma_y, double p_y);
    void set_potential(int n_slits);


};*/
/*class Matrix
{
  public:

    // Member variables
    double h, dt, T;                // Spatial and time step length and total time
    int M, N;                       // No. of points along x, y and time axes

    arma::cx_double r;              // Predefined constant, r = i∆t/h^2

    arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
    arma::cx_mat U;                 // State matrix
    arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
    arma::mat V;                    // Potential matrix

    arma::cx_cube S;                // Storing states

    // Constructor
    Matrix(double h_in, double dt_in, double T_in);

    // Takes indices (i, j) and returns the corresponding k index
    int index_k(int i, int j);

    // Set potential from file
    void set_potential(std::string filename);

    // Fill matrices A and B
    void fill_matrices();

    // Solve matrix equation Au^(n+1) = Bu^n
    void solve();

    // Set the initial state of the system
    void set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                           double sigma_y, double p_y);

};*/

/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;
        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;
        arma::cube Re;
        arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 
        double h;               
        int N;                  
        arma::cx_double r;         
        arma::sp_cx_mat A, B;          
        arma::cx_mat U;             
        arma::cx_vec u_new, u_current; 
        arma::cx_cube S;                
};*/



/*class Crank {
public:
    // Constructor
    Crank(double h, double dt, double T, double xc, double sigmax, double px,
          double yc, double sigmay, double py, double v0);

    // Function to set up the potential matrix V
    void setupPotential();

    // Function to set up the initial state matrix U0
    void setupInitialState();

    // Function to set up matrices A and B for Crank-Nicolson
    void setupMatrices();

    // Function to run the simulation loop
    void runSimulation();

    // Member variables (public for simplicity; adjust as needed)
    double h, dt, T, xc, sigmax, px, yc, sigmay, py, v0;

    // Matrices
    arma::cx_cube U;  // Assuming complex values for the wave function
    arma::mat V, A, B;
    int M = 1. / h + 1;
    int Nt = T / dt + 1;
    arma::cx_cube Re;
    arma::cx_cube Im;

    // Other helper functions if needed
};*/




/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;
        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;
        arma::cube Re;
        arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 
        double h;               
        int N;                  
        arma::cx_double r;         
        arma::sp_cx_mat A, B;          
        arma::cx_mat U;             
        arma::cx_vec u_new, u_current; 
        arma::cx_cube S;                
};*/

/*class Diff {
public:
    arma::cube storage;
    arma::cube Re;
    arma::cube Im;
    int n_timesteps;
    double T;
    double dt;
    int len;
    int M;
    double h;

    Diff(double h, double dt, double T, double xc, double yc, double px, double py, double width_x, double width_y, double potential, int n_slit = 2, double thick_x = 0.02, double center_x = 0.5, double slit_sep = 0.05, double aperture = 0.05);

    inline int index2Dto1D(int i, int j, int M);

    void generateMatrices(double r, arma::mat& A, arma::mat& B);
    void crankNicolsonStep(arma::mat& A, arma::mat& B, arma::cx_mat& u_n);
    void initializeInitialState(arma::cx_mat& u_0, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double h);
    void initializePotential(arma::mat& V, double h, double slitWidth, double wallThickness, double wallPosition, double slitSeparation, double barrierHeight);
    void run();
};*/


/*class Diff {
public:
    arma::cube storage;
    arma::cube Re;
    arma::cube Im;
    int n_timesteps;
    double T;
    double dt;
    int len;
    int M;
    Diff(double h, double dt, double T, double xc, double yc, double px, double py, double width_x, double width_y, double potential, int n_slit = 2, double thick_x = 0.02, double center_x = 0.5, double slit_sep = 0.05, double aperture = 0.05);
    inline int index2Dto1D(int i, int j, int M);
    void generateMatrices(int M, double r, arma::mat& A, arma::mat& B);
    void crankNicolsonStep(arma::mat& A, arma::mat& B, arma::vec& u_n);
    void initializeInitialState(arma::cx_mat& u_0, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double h);
    void initializePotential(arma::mat& V, double h, double slitWidth, double wallThickness, double wallPosition, double slitSeparation, double barrierHeight);
    void run();
};*/

/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;
        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;
        arma::cube Re;
        arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 



        double h;                // Spatial and time step length and total time
        int N;                       // No. of points along x, y and time axes

        arma::cx_double r;              // Predefined constant, r = i∆t/h^2

        arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
        arma::cx_mat U;                 // State matrix
        arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
                    // Potential matrix

        arma::cx_cube S;                // Storing states
        // Takes indices (i, j) and returns the corresponding k index
        
};*/












///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def///////////////////////////Two def

/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;
        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;
        arma::cube Re;
        arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 



        double h;                // Spatial and time step length and total time
        int N;                       // No. of points along x, y and time axes

        arma::cx_double r;              // Predefined constant, r = i∆t/h^2

        arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
        arma::cx_mat U;                 // State matrix
        arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
                    // Potential matrix

        arma::cx_cube S;                // Storing states
        // Takes indices (i, j) and returns the corresponding k index
        
};*/





























///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
/*class Diff
{
  public:

    // Member variables
    double h, dt, T;                // Spatial and time step length and total time
    int M, N;                       // No. of points along x, y and time axes

    arma::cx_double r;              // Predefined constant, r = i∆t/h^2

    arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
    arma::cx_mat U;                 // State matrix
    arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
    arma::mat V;                    // Potential matrix

    arma::cx_cube S;                // Storing states

    // Constructor
    Diff(double h_in, double dt_in, double T_in);

    // Takes indices (i, j) and returns the corresponding k index
    int index_k(int i, int j);

    // Set potential from file
    void set_potential(std::string filename);

    // Fill matrices A and B
    void fill_matrices();

    // Solve matrix equation Au^(n+1) = Bu^n
    void solve();

    // Set the initial state of the system
    void set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                           double sigma_y, double p_y);

};*/

/*
class Diff{
    public:
        //arma::cx_dvec u;
        arma::mat V;
        //vector<arma::sp_cx_dmat> AB;
        //arma::cx_dvec b;
        //arma::cube storage;
        double T;
        double dt;
        //int n_timesteps;
        //int len;
        int M;
        //arma::cube Re;
        //arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 
        double h;                // Spatial and time step length and total time
        int N;                       // No. of points along x, y and time axes

        arma::cx_double r;              // Predefined constant, r = i∆t/h^2

        arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
        arma::cx_mat U;                 // State matrix
        arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
                    // Potential matrix

        arma::cx_cube S;                // Storing states
        // Takes indices (i, j) and returns the corresponding k index
        int index_k(int i, int j);

        // Set potential from file
        //void set_potential(std::string filename);

        // Fill matrices A and B
        void fill_matrices();

        // Solve matrix equation Au^(n+1) = Bu^n
        //void solve();

        // Set the initial state of the system
        //void set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                            //double sigma_y, double p_y);



};*/

/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;
        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;
        arma::cube Re;
        arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 



        double h;                // Spatial and time step length and total time
        int N;                       // No. of points along x, y and time axes

        arma::cx_double r;              // Predefined constant, r = i∆t/h^2

        arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
        arma::cx_mat U;                 // State matrix
        arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
                    // Potential matrix

        arma::cx_cube S;                // Storing states
        // Takes indices (i, j) and returns the corresponding k index
        int index_k(int i, int j);

        // Set potential from file
        void set_potential(std::string filename);

        // Fill matrices A and B
        void fill_matrices();

        // Solve matrix equation Au^(n+1) = Bu^n
        void solve();

        // Set the initial state of the system
        void set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                            double sigma_y, double p_y);
};*/




/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;
        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;
        arma::cube Re;
        arma::cube Im;
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
    
        ///////////////////
        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 
};*/



/*class Diff{
    public:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;

        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;

        arma::cube Re;
        arma::cube Im;

    
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
        
        ///////////////////

        void run();
        void print_u(int i);
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len);

};*/


/*class Diff{
    private:
        arma::cx_dvec u;
        arma::mat V;
        vector<arma::sp_cx_dmat> AB;
        arma::cx_dvec b;
        arma::cube storage;

        double T;
        double dt;
        int n_timesteps;
        int len;
        int M;

        arma::cube Re;
        arma::cube Im;

    public:
        Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
        arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);
        arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);
        arma::mat probability(arma::cx_dvec &u);
        
        ///////////////////

        void run();
        void print(std::string filename);
        void print_potential(std::string filename);
        void print_u(int i);

        void save_u(std::string filename);
        
        arma::sp_cx_dmat create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i);
        arma::sp_cx_dmat create_rdiag(const std::complex<double> r, const int len);
        static int pair_to_single(const int i, const int j, const int len=0); 
        std::vector<arma::sp_cx_dmat> create_AB(arma::mat &V, const double h, const double dt, const int M);
        std::tuple<int, int> single_to_pair(const int k, const int len); 
        arma::sp_cx_dmat create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len); 



};*/


#endif