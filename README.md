# Project5 FYS3150 : Wave-Particle Duality

This repository contains c++ code constructing a class definition of the "double slit experiment", as function of an arbitrary amount of slits in a wall. This class is used to run simulations for varying amount of slits. The reopository is devided into maps `src`, `include`, `data`, `python` and `report`. Outside of these maps, in our main repository, the c++ files that run different kinds of simulations are stored. The model itself is defined in two seperate files, `diffraction.cpp` and `diffraction.hpp`, stored in maps `src` (for source files .cpp) and  `include` (for header files .hpp) dependent on their filetype. When running any simulations provided by the code in this repository, these maps have to be linked together when compiling.


----------------------
## The Class Definition of the 2D Time-Dependent Schrodinger Equation
The Iclass Diff is declared in the header file `diffraction.hpp`, and further defined in our source file `diffraction.cpp`. The purpouse of these definitions is to model the system in which a Gaussian wave packet interacts with a wall containing an arbitrary amount of slits, and how the probability density changes as time evolves in such a system. The class is modelled on the basis of a simplified environment, with x and y both being $\epsilon [0,1]$, and strict Dirichlet boundary conditions. The wall containing slits in which the wave-functions are interacting with, is modelled by implementing a potential which only lets the wave-function pass through the slits, and are reflected otherwise. 


### `diffraction.hpp`:
Contains a class `Diff` that defines and models the behaviour of wave-functions in our environment, that is also defined by functions in our class, as it evolves in time $t\epsilon [0,T]$.
The header file defines the methods:
 * `indexes(int i, int j)`, `AB_def()`, `initialize(double xc, double sigmax, double px, double yc, double sigmay, double py)`, `pot_def(int slits, int V0)` defining the environment in which simualtions are to be ran
 * `simulate()` runs simulations making use of the functions and methods defined previously

### `Idiffrection.cpp`:
Defines the methods and functions for the IsingModel class presented above: 
 * `indexes(int i, int j)` used to navigate through our 2D-matrices, implements another variable k that lets us navigate using only a single variable instead of the indices (i,j)
 * `AB_def()` defining matrices A and B needed to perform Crank-Nicolson. Needs to be called for a system defined by this class in order to run a simulation
 * `initialize(double xc, double sigmax, double px, double yc, double sigmay, double py)` initializes the system environment, by defining the initial condition of the wave-packets. This is modelled as a Gaussian wave-packet. Needs to be called for a system defined by this class in order to run a simulation
 * `pot_def(int slits, int V0)` defines the potential as a function of the amount of slits we wish the wall to contain. Needs to be called for a system defined by this class in order to run a simulation
 * `simulate()` runs simulations by making use of the Crank-Nicolson scheme



----------------------

We make use of this class when performing all of the following simualtions. All data collected by these simulations is stored as .txt files in the directory '/data'. 

### `simulation.cpp`:
Simulates and logs the behaviour of five different instances of the class Diff. Simulation1 simualtes a system where T=0.008, and there is no wall present for the wave-function to interact with. Simulation2 simualtes a system where T=0.008, and the wall the wave-fuinction interacts with contains 2 slits. The simulation3, simulation4 and simulation5 runs simulations for T=0.002 for a single, double and tripple slits. 

--------------------

## Linking and Compiling of Our Project Files:
Each simulation is ran by running the following commands in the terminal, while beinf located in the correct directory:

```sh
g++ simulation.cpp src/*.cpp -I include -o simulation -O2 -llapack -lblas -larmadillo

```
```sh
./simulation
```

----------------------

## Visualizing Data By Running Python Programs:
When ran, the simulations produce sets of data, which can be imported and visualized by running the python programs with the similar names, stored in the map 'python' in our repository: 

### `diffraction_plots.py`:
Imports the data created by `simulation.cpp`, and creates a plot visualizing the difference in deviation in probaibility from 1 for simulation1 and simulation2, as well as probability distribution along the y-axis at x=0.8 for simulation3, simulation4 and simulation5 and snapshots showcasing the probability density for simulatiion3 at t=0, t=0.001 and t=0.002.

### diffraction_animations.py`:
Imports the data created by `simulation.cpp`, and creates an animation visualizing how the probability density changes with time for the instance of each simulation, and saves the animation as a GIF. 



