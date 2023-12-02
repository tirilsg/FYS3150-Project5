import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.signal import find_peaks




#===========================================================================================================
#---------------------------- ANIMATION DEFNINTION, DEVELOPED FROM CODE EXAMPLES --------------------------#
#===========================================================================================================
def animation(filename, gifname):
    h=.005
    dt=2.5e-5
    solution = pa.cx_cube()
    solution.load(filename)
    S = np.array(solution)
    X = np.arange(0, 1 + h, h)
    Y = np.arange(0, 1 + h, h)
    x, y = np.meshgrid(X, Y)
    t = np.arange(0, 1 + dt, dt)
    t_min = t[0]
    x_min, x_max = X[0], X[-1]
    y_min, y_max = Y[0], Y[-1]
    fig = plt.figure()
    ax = plt.gca()
    img = ax.imshow(abs(S[0,:,:])**2, extent=[x_min,x_max,y_min,y_max], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(abs(S[0]))**2))
    plt.xlabel("x")
    plt.ylabel("y")
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("Probability")
    cbar.ax.tick_params()
    time_txt = plt.text(0.95, 0.95, "t = {:.3f}".format(t_min), color="white", horizontalalignment="right", verticalalignment="top")
    def animation(i):
        img.set_data(abs(S[i,:,:])**2)
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3f}".format(current_time))
        return img
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(S[:,0,0]), 2), repeat=True, blit=0)
    anim.save(gifname, writer="ffmpeg", bitrate=10000, fps=150)
    plt.show()
    
animation("../data/Simulation_1.bin", "../report/figures/Simulation1.gif")
animation("../data/Simulation_2.bin", "../report/figures/Simulation2.gif")
animation("../data/Simulation_3.bin", "../report/figures/Simulation3.gif")
animation("../data/Simulation_4.bin", "../report/figures/Simulation4.gif")
animation("../data/Simulation_5.bin", "../report/figures/Simulation5.gif")