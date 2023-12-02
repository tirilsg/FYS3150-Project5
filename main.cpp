#include "diffraction.hpp"


using namespace std;
using namespace arma;
int main(){
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T= 0.008; 
    double xc = 0.25;
    double sigmax = 0.05; 
    double px = 200; 
    double yc = 0.5;
    double sigmay = 0.05;
    double py = 0; 
    int slits = 2;
    double V0=0;
    Diff sim1 = Diff(h, dt, T);
    sim1.set_potential(slits, V0); 
    sim1.fill_matrices();
    sim1.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim1.simulate();
    sim1.S.save("data/Simulationn_1.bin");
    
    sigmay = 0.1;
    V0 = 1e10;
    Diff sim2 = Diff(h, dt, T);
    sim2.set_potential(slits, V0); 
    sim2.fill_matrices();
    sim2.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim2.simulate();
    sim2.S.save("data/Simulationn_2.bin");

    sigmay = 0.2;
    T=0.002;
    Diff sim3 = Diff(h, dt, T);
    sim3.set_potential(slits, V0); 
    sim3.fill_matrices();
    sim3.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim3.simulate();
    sim3.S.save("data/Simulationn_3.bin");

    slits = 1;
    Diff sim4 = Diff(h, dt, T);
    sim4.set_potential(slits, V0); 
    sim4.fill_matrices();
    sim4.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim4.simulate();
    sim4.S.save("data/Simulationn_4.bin");

    slits = 3;
    Diff sim5 = Diff(h, dt, T);
    sim5.set_potential(slits, V0); 
    sim5.fill_matrices();
    sim5.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim5.simulate();
    sim5.S.save("data/Simulationn_5.bin");
    return 0;
}



/*int main(){
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T= 0.008; 
    double xc = 0.25;
    double sigmax = 0.05; 
    double px = 200; 
    double yc = 0.5;
    double sigmay = 0.05;
    double py = 0; 
    int slits = 2;
    double V0=0;
    Diff sim1 = Diff(h, dt, T);
    sim1.set_potential(slits, V0); 
    sim1.fill_matrices();
    sim1.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim1.simulate();
    sim1.S.save("data/Simulation_1.bin");
    
    sigmay = 0.1;
    V0 = 1e10;
    Diff sim2 = Diff(h, dt, T);
    sim2.set_potential(slits, V0); 
    sim2.fill_matrices();
    sim2.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim2.simulate();
    sim2.S.save("data/Simulation_2.bin");

    sigmay = 0.2;
    T=0.002;
    Diff sim3 = Diff(h, dt, T);
    sim3.set_potential(slits, V0); 
    sim3.fill_matrices();
    sim3.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim3.simulate();
    sim3.S.save("data/Simulation_3.bin");

    slits = 1;
    Diff sim4 = Diff(h, dt, T);
    sim4.set_potential(slits, V0); 
    sim4.fill_matrices();
    sim4.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim4.simulate();
    sim4.S.save("data/Simulation_4.bin");

    slits = 3;
    Diff sim5 = Diff(h, dt, T);
    sim5.set_potential(slits, V0); 
    sim5.fill_matrices();
    sim5.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim5.simulate();
    sim5.S.save("data/Simulation_5.bin");
    return 0;
}*/










///////////////////////////////////////////////////////////////////////NOTES
/*int main(){
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T= 0.002; 
    double xc = 0.25;
    double sigmax = 0.05; 
    double px = 200; 
    double yc = 0.5;
    double sigmay = 0.05;
    double py = 0; 
    int slits =2;
    double V0=0;
    Diff sim1 = Diff(h, dt, T);
    sim1.set_potential(slits); 
    sim1.fill_matrices();
    sim1.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim1.solve();
    sim1.S.save("data/Simulation_1.bin");
    
    sigmay = 0.1;
    V0 = 1e10;
    Diff sim2 = Diff(h, dt, T);
    sim2.set_potential(slits); 
    sim2.fill_matrices();
    sim2.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim2.solve();
    sim2.S.save("data/Simulation_2.bin");

    sigmay = 0.2;
    T=0.002;
    Diff sim3 = Diff(h, dt, T);
    sim3.set_potential(slits); 
    sim3.fill_matrices();
    sim3.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim3.solve();
    sim3.S.save("data/Simulation_3.bin");

    slits = 1;
    Diff sim4 = Diff(h, dt, T);
    sim4.set_potential(slits); 
    sim4.fill_matrices();
    sim4.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim4.solve();
    sim4.S.save("data/Simulation_4.bin");

    slits = 3;
    Diff sim5 = Diff(h, dt, T);
    sim5.set_potential(slits); 
    sim5.fill_matrices();
    sim5.set_initial_state(xc, sigmax, px, yc, sigmay, py);
    sim5.solve();
    sim5.S.save("data/Simulation_5.bin");
    return 0;
}*/



/*int main() {

    Diff diff_instance(pass constructor parameters );

    arma::mat V = diff_instance.potential();
    arma::cx_mat u = diff_instance.wave_init();
    arma::cx_mat A, B;
    double r = ;
    diff_instance.create_AB(A, B, r);

 
    arma::cx_vec U = diff_instance.solve_system(A, B);

    // Perform other simulation steps...

    return 0;
}*/


/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.002; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 
    int n_slit = 2;
    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_1.run();
    exp_1.storage.save("data/Experimenttttt_1.bin");
    exp_1.V.save("data/Experimenttttt_1_pot.bin");

    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_2.run();
    exp_2.storage.save("data/Experimenttttt_2.bin");
    exp_2.V.save("data/Experimenttttt_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_3.run();
    exp_3.storage.save("data/Experimenttttt_3.bin");
    exp_3.Re.save("data/Experimenttttt_3_Re.bin");
    exp_3.Im.save("data/Experimenttttt_3_Im.bin");
    exp_3.V.save("data/Experimenttttt_3_pot.bin");

    n_slit = 1;
    Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_4.run();
    exp_4.storage.save("data/Experimenttttt_4.bin");
    exp_4.V.save("data/Experimenttttt_4_pot.bin");

    n_slit = 3;
    Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_5.run();
    exp_5.storage.save("data/Experimenttttt_5.bin");
    exp_5.V.save("data/Experimenttttt_5_pot.bin");
    return 0;
}*/









/*int main()
{
  int slits = 2;
  double h_in = 0.005;
  double dt_in = 2.5e-5; 
  double T_in = 0.002; 
  double x_c = 0.25; 
  double sigma_x = 0.05; 
  double p_x = 200;
  double y_c = 0.5; 
  double sigma_y = 0.05; 
  double p_y = 0; 
  sim sim = sim(h_in, dt_in, T_in, x_c, sigma_x, p_x, y_c, sigma_y, p_y, slits);
    sim.solve();
    sim.S.save("data/Experimenttt_2.bin");
    sim.V.save("data/Experimenttt_2_pot.bin");*/

  /*sim sim = sim(h_in, dt_in, T_in);
  sim.set_potential(potential);
  sim.fill_matrices();
  sim.set_initial_state(x_c, sigma_x, p_x, y_c, sigma_y, p_y);
  sim.solve();
  sim.S.save("data/Experimentt_2.bin");
  sim.V.save("data/Experimentt_2_pot.bin");*/

  //return 0;
//}
/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 

    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_1.run();
    exp_1.storage.save("data/Experimentt_1.bin");
    exp_1.V.save("data/Experimentt_1_pot.bin");

    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_2.run();
    exp_2.storage.save("data/Experimentt_2.bin");
    exp_2.V.save("data/Experimentt_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_3.run();
    exp_3.storage.save("data/Experimentt_3.bin");
    exp_3.Re.save("data/Experimentt_3_Re.bin");
    exp_3.Im.save("data/Experimentt_3_Im.bin");
    exp_3.V.save("data/Experimentt_3_pot.bin");

    int n_slit = 1;
    Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_4.run();
    exp_4.storage.save("data/Experimentt_4.bin");
    exp_4.V.save("data/Experimentt_4_pot.bin");

    n_slit = 3;
    Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_5.run();
    exp_5.storage.save("data/Experimentt_5.bin");
    exp_5.V.save("data/Experimentt_5_pot.bin");
    return 0;
}*/

/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double sigmax = 0.05; 
    double sigmay = 0.05;
    double v0 = 0; 

    Crank crank1(h, dt, T, xc, sigmax, px, yc, sigmay, py, v0);
    crank1.runSimulation();
    crank1.U.save("data/Experimentt_1.bin");
    crank1.V.save("data/Experimentt_1_pot.bin");

    v0 = 1e10;
    sigmay = 0.1; 
    Crank crank2(h, dt, T, xc, sigmax, px, yc, sigmay, py, v0);
    crank2.runSimulation();
    crank2.U.save("data/Experimentt_2.bin");
    crank2.V.save("data/Experimentt_2_pot.bin");
    return 0;
}*/
/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 

    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_1.run();
    exp_1.storage.save("data/Experimentt_1.bin");
    exp_1.V.save("data/Experimentt_1_pot.bin");

    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_2.run();
    exp_2.storage.save("data/Experimentt_2.bin");
    exp_2.V.save("data/Experimentt_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_3.run();
    exp_3.storage.save("data/Experimentt_3.bin");
    exp_3.Re.save("data/Experimentt_3_Re.bin");
    exp_3.Im.save("data/Experimentt_3_Im.bin");
    exp_3.V.save("data/Experimentt_3_pot.bin");

    int n_slit = 1;
    Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_4.run();
    exp_4.storage.save("data/Experimentt_4.bin");
    exp_4.V.save("data/Experimentt_4_pot.bin");

    n_slit = 3;
    Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_5.run();
    exp_5.storage.save("data/Experimentt_5.bin");
    exp_5.V.save("data/Experimentt_5_pot.bin");
    return 0;
}*/

/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 

    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_1.run();
    exp_1.storage.save("data/Experimentt_1.bin");
    exp_1.V.save("data/Experimentt_1_pot.bin");

    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_2.run();
    exp_2.storage.save("data/Experimentt_2.bin");
    exp_2.V.save("data/Experimentt_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_3.run();
    exp_3.storage.save("data/Experimentt_3.bin");
    exp_3.Re.save("data/Experimentt_3_Re.bin");
    exp_3.Im.save("data/Experimentt_3_Im.bin");
    exp_3.V.save("data/Experimentt_3_pot.bin");

    int n_slit = 1;
    Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_4.run();
    exp_4.storage.save("data/Experimentt_4.bin");
    exp_4.V.save("data/Experimentt_4_pot.bin");

    n_slit = 3;
    Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_5.run();
    exp_5.storage.save("data/Experimentt_5.bin");
    exp_5.V.save("data/Experimentt_5_pot.bin");
    return 0;
}*/

/*int main() {
    Diff myDiff(1.0, 1.0, 100.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 1.0);  // Adjust parameters accordingly
    myDiff.run();
    myDiff.storage.save("result_cube.h5", arma::hdf5_binary);

    return 0;
}*/

/*int main() {
    double h = 0.005;
    double dt = 2.5e-5;
    double T = 0.008;
    double xc = 0.25;
    double yc = 0.5;
    double px = 200;
    double py = 0;
    double width_x = 0.05;
    double width_y = 0.1; // corrected as per your comment
    double potential = 1e10;

    Diff diffraction(1.0, 1.0, 100, 0.1, 0.2, 0.3, 0.4);
    diffraction.run();

    // Save the result
    myDiff.storage.save("result_cube.h5", arma::hdf5_binary);
    return 0;
}*/

/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 

    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    //Diff exp_1 = Diff(h, dt, T);
    //exp_1.fill_matrices();
    cout << exp_1.V << endl ;
    cout << exp_1.A << endl ;

    //Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    //exp_1.run();
    //exp_1.storage.save("data/Experiment_1.bin");

    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    //exp_2.fill_matrices();
    cout << exp_2.V << endl ;
    cout << exp_2.A << endl ;
    exp_2.run();
    exp_2.storage.save("data/Experiment_2.bin");
    exp_2.V.save("data/Experiment_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    //Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    //exp_3.fill_matrices();
    //cout << exp_3.V << endl ;
    //exp_3.run();
    //exp_3.storage.save("data/Experiment_3.bin");
    //exp_3.Re.save("data/Experiment_3_Re.bin");
    //exp_3.Im.save("data/Experiment_3_Im.bin");

    int n_slit = 1;
    //Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    //exp_4.fill_matrices();
    //cout << exp_4.V << endl ;
    //exp_4.run();
    //exp_4.storage.save("data/Experiment_4.bin");

    n_slit = 3;
    //Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    //exp_5.fill_matrices();
    //cout << exp_5.V << endl ;
    //exp_5.run();
    //exp_5.storage.save("data/Experiment_5.bin");
  
    return 0;
}*/

/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 

    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_1.fill_matrices();
    exp_1.set_initial_state(xc, width_x, px, yc, width_y, py);
    exp_1.solve();
    exp_1.S.save("data/Experiment_1.bin");



    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_2.fill_matrices();
    exp_2.set_initial_state(xc, width_x, px, yc, width_y, py);
    exp_2.solve();
    exp_2.S.save("data/Experiment_2.bin");
    exp_2.V.save("data/Experiment_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_3.fill_matrices();
    exp_3.set_initial_state(xc, width_x, px, yc, width_y, py);
    exp_3.solve();
    exp_3.S.save("data/Experiment_3.bin");
    exp_3.Re.save("data/Experiment_3_Re.bin");
    exp_3.Im.save("data/Experiment_3_Im.bin");

    int n_slit = 1;
    Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_4.fill_matrices();
    exp_4.set_initial_state(xc, width_x, px, yc, width_y, py);
    exp_4.solve();
    exp_4.S.save("data/Experiment_4.bin");

    n_slit = 3;
    Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_5.fill_matrices();
    exp_5.set_initial_state(xc, width_x, px, yc, width_y, py);
    exp_5.solve();
    exp_5.S.save("data/Experiment_5.bin");
  
    return 0;
}*/



/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 

    //Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    //exp_1.fill_matrices();
    //exp_1.set_initial_state(xc, width_x, px, yc, width_y, py);
    //exp_1.solve();
    //exp_1.S.save("data/Experiment_1.bin");

    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_1.run();
    exp_1.storage.save("data/Experiment_1.bin");

    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_2.run();
    exp_2.storage.save("data/Experiment_2.bin");
    exp_2.V.save("data/Experiment_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_3.run();
    exp_3.storage.save("data/Experiment_3.bin");
    exp_3.Re.save("data/Experiment_3_Re.bin");
    exp_3.Im.save("data/Experiment_3_Im.bin");

    int n_slit = 1;
    Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_4.run();
    exp_4.storage.save("data/Experiment_4.bin");

    n_slit = 3;
    Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_5.run();
    exp_5.storage.save("data/Experiment_5.bin");
  
    return 0;
}*/



/*int main() {
    double h = 0.005; 
    double dt = 2.5e-5; 
    double T = 0.008; 
    double xc = 0.25; 
    double yc = 0.5; 
    double px = 200; 
    double py = 0; 
    double width_x = 0.05; 
    double width_y = 0.05;
    double potential = 0; 

    Diff exp_1 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_1.run();
    exp_1.storage.save("data/Experiment_1.bin");

    potential = 1e10;
    width_y = 0.1; 
    Diff exp_2 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_2.run();
    exp_2.storage.save("data/Experiment_2.bin");
    exp_2.V.save("data/Experiment_2_pot.bin");

    width_y = 0.2; 
    T = 0.002; 
    Diff exp_3 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential);
    exp_3.run();
    exp_3.storage.save("data/Experiment_3.bin");
    exp_3.Re.save("data/Experiment_3_Re.bin");
    exp_3.Im.save("data/Experiment_3_Im.bin");

    int n_slit = 1;
    Diff exp_4 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_4.run();
    exp_4.storage.save("data/Experiment_4.bin");

    n_slit = 3;
    Diff exp_5 = Diff(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    exp_5.run();
    exp_5.storage.save("data/Experiment_5.bin");
  
    return 0;
}*/