#python program that imports data files collected by compiling and running simulation defined in simulation.cpp
#collects data of matrices gathered in a cube, visualizing how the wave-function in the 2d xy-space evolves with time, in intervals of time T=0.008 and 0.002
#uses this data to vislualize the evolution in time, in the shape of contour in an animated GIF. These GIFs are probability plots
#each animation is saved in directory "../report/figures/", which is consistent with compilation of latex file to pdf for the report
#this program is based on the code example provided 
#we start by importing needed modules
import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#===========================================================================================================
#---------------------------- ANIMATION DEFNINTION, DEVELOPED FROM CODE EXAMPLES --------------------------#
#===========================================================================================================
def animation(filename, gifname): #defines a function that animates the dataset, as a function of direction to data and direction where we wish to save the animation
    h=0.005
    dt=2.5e-5
    solution = pa.cx_cube()
    solution.load(filename) #loads the data
    S = np.array(solution) #creates an array containing all the data
    X = np.arange(0, 1 + h, h) #defines the xy-space
    Y = np.arange(0, 1 + h, h) #defines the xy-space
    x, y = np.meshgrid(X, Y)
    t = np.arange(0, 1 + dt, dt) #defines the t-space
    #defines a couple "variables"
    t_min = t[0]
    x_min, x_max = X[0], X[-1]
    y_min, y_max = Y[0], Y[-1]
    fig = plt.figure() #initializes a figure, and susequently defines the animation as in the codeexample
    ax = plt.gca()
    #color scale normalization 
    img = ax.imshow(abs(S[0,:,:])**2, extent=[x_min,x_max,y_min,y_max], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(abs(S[0]))**2))
    plt.xlabel("x")
    plt.ylabel("y")
    #colorbar implementation
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("Probability")
    cbar.ax.tick_params()
    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3f}".format(t_min), color="white", horizontalalignment="right", verticalalignment="top")
    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        img.set_data(abs(S[i,:,:])**2)
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3f}".format(current_time))
        return img
    #creates animation, and both shows and saves it in desired directory
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(S[:,0,0]), 2), repeat=True, blit=0)
    anim.save(gifname, writer="ffmpeg", bitrate=10000, fps=150)
    plt.show()

#calls the function for all our data sets we wish to create animations for
animation("../data/Simulation_1.bin", "../report/figures/Simulation1.gif")
animation("../data/Simulation_2.bin", "../report/figures/Simulation2.gif")
animation("../data/Simulation_3.bin", "../report/figures/Simulation3.gif")
animation("../data/Simulation_4.bin", "../report/figures/Simulation4.gif")
animation("../data/Simulation_5.bin", "../report/figures/Simulation5.gif")