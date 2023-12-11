#python program that imports data files collected by compiling and running simulation defined in simulation.cpp
#collects data of matrices gathered in a cube, visualizing how the wave-function in the 2d xy-space evolves with time, in intervals of time T=0.008 and 0.002
#uses this data to create different sets of plots, including deviation, probability density plots
#each plot is saved in directory "../report/figures/", which is consistent with compilation of latex file to pdf for the report
#this program is based on the code example provided 
#we start by importing needed modules
import matplotlib
import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use("ggplot") #style of plots
cmap_magma = plt.cm.get_cmap("magma") #colors we want to use to create plots, the same color scheme as the corresponding animations
colors = [cmap_magma(i/14) for i in range(15)]

def calculate_errors(paths, names): #definition of a function calculating error for the simulations ran in the time interval T=0.008
    T=0.008
    all_errors = [] #list that is goint to contain errors
    for path, name in zip(paths, names): #loops through the data and calculates errors, and appends to the list
        solution = pa.cx_cube()
        solution.load(path) #loads data
        S = np.array(solution)
        t = np.linspace(0, T, len(S))
        deviation = np.real(np.sum(np.conj(S) * S, axis=(1, 2)))#calculates probability
        abs_err = np.abs(1 - deviation) #calculates the deviation from 1
        avrg_err = np.mean(abs_err) #expresses the average error for the entire time interval for the data set path
        all_errors.append((t, abs_err, name, avrg_err)) #appends errors to the list
    return all_errors #returns list containing errors for all paths in our path list

paths = ["../data/Simulation_1.bin", "../data/Simulation_2.bin"] #paths where the simulation data is stored
names = ["No Slits", "Double Slits"] #what data each path contains
data = calculate_errors(paths, names) #calculates errors 
i=4 #helps us decide what colors to use 
for t, abs_err, name, avrg_err in data: #plots the deviation from 1 for each data set
    plt.plot(t, abs_err, color=colors[i], linewidth=2.5, label=f"{name}:\nAverage error: {avrg_err:.3e}") #adds the average error in the legend of the plots
    i+=6 #helps us decide what colors to use 
plt.xlabel("Time [s]")
plt.ylabel("Deviation from 1")
plt.tight_layout()
plt.legend()
plt.savefig("../report/figures/deviation.pdf") #saves figure
plt.show()



# ===========================================================================================================
# -------------------------- PLOT OF POTENTIAL MEASURED BY DETECTOR SCREEN AT X=0.8 ------------------------
# ===========================================================================================================
t_index = int(0.002 / 2.5e-5) 
x_index = int(0.8 / 0.005)
def plot_probability(x_values, y_values, label_text, line_color): #defines the plotting as a function of x and y, as well as label and the desired color of the plot
    plt.plot(x_values, y_values, label=label_text, color=line_color, linewidth=2.5) #plots the functions
    peaks, _ = find_peaks(y_values, height=0.005) #finds the peaks using function imported as module
    for i, peak in enumerate(peaks): #plots the peaks of the functions with text describing its coordinates
        plt.scatter(x_values[peak], y_values[peak], color=line_color, marker="x", linewidth=2.5)
        plt.text(x_values[peak], y_values[peak], f"({x_values[peak]:.3f}, {y_values[peak]:.3f})", color=line_color)
        
plt.figure(figsize=(8, 6)) #creates figure
data = pa.cx_cube()
data.load("../data/Simulation_3.bin") #loads data
S = np.array(data) ; P = np.conj(S) * S
probability = P[t_index, :, x_index] / np.sum(P[t_index, :, x_index]) #probability function definition from data
y_values = np.linspace(0, 1, len(S[0])) # y definition at x=0.8
plot_probability(y_values, np.real(probability), "Double Slits, Simulation 3", colors[4]) #creates plots

data.load("../data/Simulation_4.bin")
S = np.array(data) ; P = np.conj(S) * S
probability = P[t_index, :, x_index] / np.sum(P[t_index, :, x_index])
y_values = np.linspace(0, 1, len(S[0]))
plot_probability(y_values, np.real(probability), "Single Slit, Simulation 4", colors[10])

data.load("../data/Simulation_5.bin")
S = np.array(data) ; P = np.conj(S) * S
probability = P[t_index, :, x_index] / np.sum(P[t_index, :, x_index])
y_values = np.linspace(0, 1, len(S[0]))
plot_probability(y_values, np.real(probability), "Triple Slits, Simulation 5", colors[-3])
plt.xlabel("y")
plt.ylabel("Probability")
plt.legend()
plt.savefig("../report/figures/wall.pdf") #saves figure
plt.show()



#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#----------------------------------- BOTH IMAGINARY AND REAL DIMENSIONS -----------------------------------#
#===========================================================================================================
dt = 2.5e-5
solution = pa.cx_cube()
solution.load("../data/Simulation_3.bin") #loads data 
t_list = np.array([0.0, 0.001, 0.002])
indexes=(t_list/dt)

def plots(func, name): #function creating the contour plot
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for j in range(3): #defines three plots
        axes[j].text(0.9, 0.9, f"t = {t_list[j]:.3f}", color="white", horizontalalignment="right", verticalalignment="top") #adds text telling what time the snapshot is visualizing
        im = axes[j].imshow(func[int(indexes[j]), :, :], extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(func[0]))) #creates probability density plot
    divider = make_axes_locatable(axes[2])
    cax = divider.append_axes("right", size="5%", pad=0.05) #makes sure all three plots have the same size, not affected by the affed colorbar
    cbar = fig.colorbar(im, cax=cax, label="Probability") #adds colorbar with label
    axes[1].set_xlabel("x-dimension")
    axes[0].set_ylabel("y-dimension")
    axes[0].grid() #removes grid 
    axes[1].grid() #removes grid 
    axes[2].grid() #removes grid 
    plt.savefig(f"{name}") #saves figure

S= abs(np.array(solution)) 
plots(S, "../report/figures/plots_both.pdf") #creates plot for both real and imaginary dimensions
plt.show()

#creating 3D plot for t=0.002 using the same data, for double slits
fig = plt.figure(figsize=(8, 4)) 
ax = fig.add_subplot(111, projection="3d")
X, Y = np.meshgrid(np.linspace(0, 1, len(S[int(0.002 / dt), :, :])), np.linspace(0, 1, len(S[int(0.002 / dt), :, :])))
grr = ax.plot_surface(X, Y, S[int(0.002 / dt), :, :], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S[0])))
ax.set_xlabel("x-dimension")
ax.set_ylabel("y-dimension")
ax.grid()
cbar = fig.colorbar(grr, ax=ax, label="Probability")
plt.savefig(f"../report/figures/plot_both_3d.pdf", bbox_inches="tight") #saves figure
plt.show()

#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#--------------------------------------------- REAL DIMENSIONS --------------------------------------------#
#===========================================================================================================
S = np.real(np.array(solution))
plots(S, "../report/figures/plots_real.pdf")
plt.show()

#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#------------------------------------------- IMAGINARY DIMENSIONS -----------------------------------------#
#===========================================================================================================
S = np.imag(np.array(solution))
plots(S, "../report/figures/plots_imag.pdf")
plt.show()