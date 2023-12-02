import matplotlib
import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use("ggplot")
cmap_magma = plt.cm.get_cmap('magma')
colors = [cmap_magma(i/14) for i in range(15)]

def calculate_deviation(S_matrix):
    N_samples = len(S_matrix)
    deviations = np.zeros(N_samples)
    for n in range(N_samples):
        S_n = S_matrix[n, :, :]
        deviation = np.real(np.sum(np.conj(S_n) * S_n))
        deviations[n] = deviation
    return deviations

def calculate_errors(deviation_values):
    return abs(1 - deviation_values)

def plot_probability(x_values, y_values, label_text, line_color):
    plt.plot(x_values, y_values, label=label_text, color=line_color, linewidth=2.5)
    peaks, _ = find_peaks(y_values, height=0.005)
    for i, peak in enumerate(peaks):
        plt.scatter(x_values[peak], y_values[peak], color=line_color, marker='x', linewidth=2.5)
        plt.text(x_values[peak], y_values[peak], f'({x_values[peak]:.3f}, {y_values[peak]:.3f})', color=line_color)

# ===========================================================================================================
# ------------------------- PLOT VISUALIZING DEVIATION FOR NO BARRIERS VS. BARRIERS ------------------------
# ===========================================================================================================
T_time = 0.008
solution_no_barrier = pa.cx_cube()
solution_barrier = pa.cx_cube()
solution_no_barrier.load("../data/Simulation_1.bin")
solution_barrier.load("../data/Simulation_2.bin")

S_no_barrier = np.array(solution_no_barrier)
S_barrier = np.array(solution_barrier)
time_no_barrier = np.linspace(0, T_time, len(S_no_barrier))
time_barrier = np.linspace(0, T_time, len(S_barrier))
deviations_no_barrier = calculate_deviation(S_no_barrier)
deviations_barrier = calculate_deviation(S_barrier)
absolute_errors_no_barrier = calculate_errors(deviations_no_barrier)
absolute_errors_barrier = calculate_errors(deviations_barrier)

plt.figure(figsize=(8,7))
plt.plot(time_no_barrier, absolute_errors_no_barrier, color=colors[4], linewidth=2.5, label=f"No Slits:\nInitial error: {absolute_errors_no_barrier[0]:.3e}\nFinal error: {absolute_errors_no_barrier[-1]:.3e}")
plt.plot(time_barrier, absolute_errors_barrier, color=colors[10], linewidth=2.5, label=f"Double Slits:\nInitial error: {absolute_errors_barrier[0]:.3e}\nFinal error: {absolute_errors_barrier[-1]:.3e}")
plt.xlabel("Time [s]:")
plt.ylabel("Relative Error:")
plt.tight_layout()
plt.legend()
plt.savefig("../report/figures/deviation.pdf")
plt.show()
"""print("Time [s]:")
print(time_no_barrier)
print("Relative Error no Barrier:")
print(S_no_barrier)
print("Relative Error Barrier:")
print(S_barrier)"""
# ===========================================================================================================
# -------------------------- PLOT OF POTENTIAL MEASURED BY DETECTOR SCREEN AT X=0.8 ------------------------
# ===========================================================================================================
plt.figure(figsize=(8, 6))
data_detector_screen = pa.cx_cube()
data_detector_screen.load('../data/Simulation_3.bin')
S_detector_screen = np.array(data_detector_screen)
t_index, x_index = int(0.002 / 2.5e-5), int(0.8 / 0.005)
P_values = np.conj(S_detector_screen) * S_detector_screen
probability_values = P_values[t_index, :, x_index] / np.sum(P_values[t_index, :, x_index])
y_values = np.linspace(0, 1, len(S_detector_screen[0]))
plot_probability(y_values, np.real(probability_values), "Double Slits, Simulation 3", colors[4])

data_detector_screen.load('../data/Simulation_4.bin')
S_detector_screen = np.array(data_detector_screen)
P_values = np.conj(S_detector_screen) * S_detector_screen
probability_values = P_values[t_index, :, x_index] / np.sum(P_values[t_index, :, x_index])
plot_probability(y_values, np.real(probability_values), "Single Slit, Simulation 4", colors[10])

data_detector_screen.load('../data/Simulation_5.bin')
S_detector_screen = np.array(data_detector_screen)
P_values = np.conj(S_detector_screen) * S_detector_screen
probability_values = P_values[t_index, :, x_index] / np.sum(P_values[t_index, :, x_index])
plot_probability(y_values, np.real(probability_values), "Triple Slits, Simulation 5", colors[-3])
plt.xlabel('y')
plt.ylabel('Probability')
plt.legend()
plt.savefig("../report/figures/wall.pdf")
plt.show()




#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#----------------------------------- BOTH IMAGINARY AND REAL DIMENSIONS -----------------------------------#
#===========================================================================================================
dt = 2.5e-5
solution = pa.cx_cube()
solution.load('../data/Simulation_3.bin')
S = np.array(solution)
S_cx = abs(S)**2
S_real = np.real(S)
S_imag = np.imag(S)
S_list = [S_cx, S_real, S_imag]
x_ticks = np.linspace(0.2, 0.8, 4)
x_labels = ['0.2', '0.4', '0.6', '0.8']
t_list = np.array([0.0, 0.001, 0.002])
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_cx[idx, :, :]
    axes[j].text(0.95, 0.95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_cx[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
fig.suptitle("Both dimensions")
plt.savefig(f"../report/figures/plots_both.pdf")
plt.show()

fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111, projection='3d')
t_idx = int(0.002 / dt)
S_n_idx = S_cx[t_idx, :, :]
X, Y = np.meshgrid(np.linspace(0, 1, len(S_n_idx)), np.linspace(0, 1, len(S_n_idx[0])))
ax.plot_surface(X, Y, S_n_idx, cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_cx[0])))
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels)
ax.set_xlabel("x-dimension")
ax.set_ylabel("y-dimension")
ax.set_zlabel("Probability")
ax.set_title(f't = {t_list[2]:.3f}')
plt.savefig(f"../report/figures/plot_both_3d.pdf", bbox_inches='tight')
plt.show()

#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#--------------------------------------------- REAL DIMENSIONS --------------------------------------------#
#===========================================================================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_real[idx, :, :]
    axes[j].text(.95, .95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_real[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
fig.suptitle("Real dimensions")
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
plt.savefig(f"../report/figures/plots_real.pdf")
plt.show()

fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111, projection='3d')
t_idx = int(0.002 / dt)
S_n_idx = S_cx[t_idx, :, :]
X, Y = np.meshgrid(np.linspace(0, 1, len(S_n_idx)), np.linspace(0, 1, len(S_n_idx[0])))
ax.plot_surface(X, Y, S_n_idx, cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_cx[0])))
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels)
ax.set_xlabel("x-dimension")
ax.set_ylabel("y-dimension")
ax.set_zlabel("Probability")
ax.set_title(f't = {t_list[2]:.3f}')
plt.savefig(f"../report/figures/plot_real_3d.pdf", bbox_inches='tight')
plt.show()

#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#------------------------------------------- IMAGINARY DIMENSIONS -----------------------------------------#
#===========================================================================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_imag[idx, :, :]
    axes[j].text(.95, .95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_imag[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
fig.suptitle("Imaginary dimension")
plt.savefig(f"../report/figures/plots_imag.pdf")
plt.show()

fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111, projection='3d')
t_idx = int(0.002 / dt)
S_n_idx = S_cx[t_idx, :, :]
X, Y = np.meshgrid(np.linspace(0, 1, len(S_n_idx)), np.linspace(0, 1, len(S_n_idx[0])))
ax.plot_surface(X, Y, S_n_idx, cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_cx[0])))
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels)
ax.set_xlabel("x-dimension")
ax.set_ylabel("y-dimension")
ax.set_zlabel("Probability")
ax.set_title(f't = {t_list[2]:.3f}')
plt.savefig(f"../report/figures/plot_imag_3d.pdf", bbox_inches='tight')

plt.show()




"""

#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#----------------------------------- BOTH IMAGINARY AND REAL DIMENSIONS -----------------------------------#
#===========================================================================================================
dt = 2.5e-5
solution = pa.cx_cube()
solution.load('../data/Simulation_3.bin')
S = np.array(solution)
S_cx = abs(S)**2
S_real = np.real(S)
S_imag = np.imag(S)
S_list = [S_cx, S_real, S_imag]
x_ticks = np.linspace(0.2, 0.8, 4)
x_labels = ['0.2', '0.4', '0.6', '0.8']
t_list = np.array([0.0, 0.001, 0.002])
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_cx[idx, :, :]
    axes[j].text(0.95, 0.95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_cx[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
fig.suptitle("Both dimensions")
plt.savefig(f"../report/figures/plots_both.pdf")
plt.show()

#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#--------------------------------------------- REAL DIMENSIONS --------------------------------------------#
#===========================================================================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_real[idx, :, :]
    axes[j].text(.95, .95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_real[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
fig.suptitle("Real dimensions")
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
plt.savefig(f"../report/figures/plots_real.pdf")
plt.show()


#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#------------------------------------------- IMAGINARY DIMENSIONS -----------------------------------------#
#===========================================================================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_imag[idx, :, :]
    axes[j].text(.95, .95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_imag[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
fig.suptitle("Imaginary dimension")
plt.savefig(f"../report/figures/plots_imag.pdf")
plt.show()
"""
"""import matplotlib
import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use("ggplot")
cmap_magma = plt.cm.get_cmap('magma')
colors = [cmap_magma(i/14) for i in range(15)]

def calculate_deviation(S):
    N = len(S)
    devs = np.zeros(N)
    for n in range(N):
        S_n = S[n, :, :]
        dev = np.real(np.sum(np.conj(S_n) * S_n))
        devs[n] = dev
    return devs

def calculate_errors(deviation):
    return abs(1 - deviation)

def plot_probability(x, y, label, color):
    plt.plot(x, y, label=label, color=color, linewidth=2.5)
    peaks, _ = find_peaks(y, height=0.005)
    for i, peak in enumerate(peaks):
        plt.scatter(x[peak], y[peak], color=color, marker='x', linewidth=2.5)
        plt.text(x[peak], y[peak], f'({x[peak]:.3f}, {y[peak]:.3f})', color=color)

#===========================================================================================================
#------------------------- PLOT VISUALIZING DEVIATION FOR NO BARRIERS VS. BARRIERS ------------------------#
#===========================================================================================================
T = 0.008
solution1 = pa.cx_cube()
solution2 = pa.cx_cube()
solution1.load("../data/Simulation_1.bin")
solution2.load("../data/Simulation_2.bin")

S1 = np.array(solution1)
S2 = np.array(solution2)
t1 = np.linspace(0, T, len(S1))
t2 = np.linspace(0, T, len(S2))
devs1 = calculate_deviation(S1)
devs2 = calculate_deviation(S2)
abs_errs1 = calculate_errors(devs1)
abs_errs2 = calculate_errors(devs2)
plt.figure(figsize=(8, 4))

plt.plot(t1, abs_errs1, color=colors[4],linewidth=2.5, label=f'Init. err.: {abs_errs1[0]:.3e}\nFinal err.: {abs_errs1[-1]:.3e}')
plt.plot(t2, abs_errs2, color=colors[10], linewidth=2.5, label=f'Init. err.: {abs_errs2[0]:.3e}\nFinal err.: {abs_errs2[-1]:.3e}')
plt.xlabel('Time [s]')
plt.ylabel('Rel. error')
plt.tight_layout()
plt.legend()
#plt.grid()
#plt.savefig("../report/figures/deviation.pdf")
plt.show()





#===========================================================================================================
#-------------------------- PLOT OF POTENTIAL MEASURED BY DETECTOR SCREEN AT X=0.8 ------------------------#
#===========================================================================================================
plt.figure(figsize=(8, 6))
data = pa.cx_cube()
data.load('../data/Simulation_3.bin')
S = np.array(data)
t_idx, x_idx = int(0.002 / 2.5e-5), int(0.8 / 0.005)
P = np.conj(S) * S
p = P[t_idx, :, x_idx] / np.sum(P[t_idx, :, x_idx])
y = np.linspace(0, 1, len(S[0]))
plot_probability(y, np.real(p), 'Simulation 3', colors[4])

data.load('../data/Simulation_4.bin')
S = np.array(data)
P = np.conj(S) * S
p = P[t_idx, :, x_idx] / np.sum(P[t_idx, :, x_idx])
plot_probability(y, np.real(p), 'Simulation 4', colors[10])

data.load('../data/Simulation_5.bin')
S = np.array(data)
P = np.conj(S) * S
p = P[t_idx, :, x_idx] / np.sum(P[t_idx, :, x_idx])
plot_probability(y, np.real(p), 'Simulation 5', colors[-3])

plt.xlabel('y')
plt.ylabel('Probability')
plt.legend()
#plt.grid()
#plt.savefig("../report/figures/wall.pdf")
plt.show()


#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#----------------------------------- BOTH IMAGINARY AND REAL DIMENSIONS -----------------------------------#
#===========================================================================================================
dt = 2.5e-5
solution = pa.cx_cube()
solution.load('../data/Simulation_3.bin')
S = np.array(solution)
S_cx = abs(S)**2
S_real = np.real(S)
S_imag = np.imag(S)
S_list = [S_cx, S_real, S_imag]
x_ticks = np.linspace(0.2, 0.8, 4)
x_labels = ['0.2', '0.4', '0.6', '0.8']
t_list = np.array([0.0, 0.001, 0.002])
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_cx[idx, :, :]
    axes[j].text(0.95, 0.95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_cx[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
fig.suptitle("Both dimensions")
#plt.savefig(f"../report/figures/plots_both.pdf")
plt.show()

#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#--------------------------------------------- REAL DIMENSIONS --------------------------------------------#
#===========================================================================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_real[idx, :, :]
    axes[j].text(.95, .95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_real[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
fig.suptitle("Real dimensions")
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
#plt.savefig(f"../report/figures/plots_real.pdf")
plt.show()


#===========================================================================================================
#------------------------- EVOLUTION OF THE 2D PROBABILITY FUNCTION, DOUBLE SLITS -------------------------#
#------------------------------------------- IMAGINARY DIMENSIONS -----------------------------------------#
#===========================================================================================================
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for j in range(3):
    idx = int(t_list[j] / dt)
    S_n_idx = S_imag[idx, :, :]
    axes[j].text(.95, .95, f't = {t_list[j]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
    im = axes[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap="magma", norm=matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_imag[0])))
    axes[j].set_xticks(x_ticks)
    axes[j].set_xticklabels(x_labels)
divider = make_axes_locatable(axes[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=8)  
axes[1].set_xlabel("x-dimension")
axes[0].set_ylabel("y-dimension")
axes[0].grid()
axes[1].grid()
axes[2].grid()
fig.suptitle("Imaginary dimension")
#plt.savefig(f"../report/figures/plots_imag.pdf")
plt.show()"""