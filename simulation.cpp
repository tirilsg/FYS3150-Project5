#include "diffraction.hpp"

using namespace std;
using namespace arma;
int main(){
    // definition of variables needed for our simulations, defining both the initial condition of wave-function as well as the time interval 
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T= 0.008; // time interval of 0.008 seconds
    double xc = 0.25;
    double sigmax = 0.05; 
    double px = 200; 
    double yc = 0.5;
    double sigmay = 0.05;
    double py = 0; 

    // running simulations 
    int slits = 2;
    double V0=0; //V0 is set to 0
    Diff sim1 = Diff(h, dt, T); // defines a class-environment
    sim1.pot_def(slits, V0);  // defines potential 
    sim1.AB_def(); // defines A and B matrices
    sim1.initialize(xc, sigmax, px, yc, sigmay, py); // initializes system
    sim1.simulate(); // runs simulations
    sim1.S.save("data/Simulation_1.bin"); // saves the data for the wave-function simulation in a map called data

    sigmay = 0.1;
    V0 = 1e10; //V0 is set to 1e10
    Diff sim2 = Diff(h, dt, T);
    sim2.pot_def(slits, V0); 
    sim2.AB_def();
    sim2.initialize(xc, sigmax, px, yc, sigmay, py);
    sim2.simulate();
    sim2.S.save("data/Simulation_2.bin");

    sigmay = 0.2; // update sigmay 
    T=0.002; // time onterval of 0.002 seconds
    Diff sim3 = Diff(h, dt, T);
    sim3.pot_def(slits, V0); 
    sim3.AB_def();
    sim3.initialize(xc, sigmax, px, yc, sigmay, py);
    sim3.simulate();
    sim3.S.save("data/Simulation_3.bin");

    slits = 1;
    Diff sim4 = Diff(h, dt, T);
    sim4.pot_def(slits, V0); 
    sim4.AB_def();
    sim4.initialize(xc, sigmax, px, yc, sigmay, py);
    sim4.simulate();
    sim4.S.save("data/Simulation_4.bin");

    slits = 3;
    Diff sim5 = Diff(h, dt, T);
    sim5.pot_def(slits, V0); 
    sim5.AB_def();
    sim5.initialize(xc, sigmax, px, yc, sigmay, py);
    sim5.simulate();
    sim5.S.save("data/Simulation_5.bin");
    return 0;
}
