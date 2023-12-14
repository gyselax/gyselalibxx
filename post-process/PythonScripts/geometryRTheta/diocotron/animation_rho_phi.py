"""
Save and plot the animation of the density and the electrical potential of the diocotron simulation.
"""

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gysdata import DiskStore


# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent
default_output_dir = current_folder.joinpath("../../../../build/simulations/geometryRTheta/diocotron/output/").resolve()

parser = argparse.ArgumentParser(description="Plot and save the density and the electrical potential solutions of a given output folder.")

parser.add_argument('--folder', metavar='folder', type=str,
                    default = default_output_dir,
                    help=f"Path to the output folder. By default folder = '{default_output_dir}'")

parser.add_argument('-name', metavar='name', type=str, nargs= '?',
                    default = "diocotron",
                    help="Name of the saved animation. By default name = 'diocotron'.")

parser.add_argument('--perturb', action='store_true', help='Plot and save only the perturbation (f - f_{eq}).')

args = parser.parse_args()

folder = Path(args.folder)
name_animation = args.name
if_perturb = args.perturb

# ---------------------------------------------------------------------------

path_data_structure = Path('data_structure_diocotron.yaml')
ds = DiskStore(folder, data_structure=path_data_structure)


# Get initial data
dt = float(ds["delta_t"])
T = float(ds["final_T"])
iter_nb = int( T / dt ) +1
tstep_diag = int(ds["time_step_diag"])

Nr = int(ds["r_size"])
Nt = int(ds["p_size"])

omega_Im = float(ds["slope"])
jacobian = np.array(ds["jacobian"])


# Get initial data
x_grid = np.array(ds['x_coords']).ravel()
y_grid = np.array(ds['y_coords']).ravel()



# --- equiibrium data
if if_perturb:
    Rho_eq = np.array(ds['density_eq'])
    Phi_eq = np.array(ds['electrical_potential_eq'])
else:
    Rho_eq = np.zeros((Nr+3, Nt))
    Phi_eq = np.zeros((Nr+3, Nt))


# --- function data
Rho = np.array(ds['density'])
Phi = np.array(ds['electrical_potential'])

zarray_rho = np.array([[x_grid, y_grid, (Rho[i] - Rho_eq).ravel() ] for i in range(iter_nb // tstep_diag)])
zarray_phi = np.array([[x_grid, y_grid, (Phi[i] - Phi_eq).ravel() ] for i in range(iter_nb // tstep_diag)])


zarray = np.array([zarray_rho, zarray_phi])


# --- animation
fig = plt.figure(figsize=(16,8))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

elev, azim, roll = 90, -90, 0
ax1.view_init(elev, azim, roll)
ax2.view_init(elev, azim, roll)


my_cmap = plt.get_cmap('seismic')#('Spectral')#('plasma')#('jet')#('gnuplot')
def update_plot_function(figure, axis1, axis2):
    """
    Define a update_plot function which update the plots
    in the FuncAnimation function.

    Parameters
    ----------
    figure: matplotlib.pyplot.figure
        The figure object of the plots.
    axis1: matplotlib.pyplot.axis
        The first subplot object.
    axis2: matplotlib.pyplot.axis
        The second subplot object.

    Returns
    -------
        update_plot, the function which update the plots
        for the FuncAnimation function.
    """
    def update_plot(frame_number, zarray, plot1, colorbar1, plot2, colorbar2):
        vmin_rho = min(min(l) for l in zarray[0,:,2])
        vmax_rho = max(max(l) for l in zarray[0,:,2])

        plot1[0].remove()
        plot1[0] = axis1.plot_trisurf(zarray[0,frame_number,0], zarray[0,frame_number,1], zarray[0,frame_number,2],
                                  cmap = my_cmap, linewidth=0, antialiased=False,
                                  vmin = vmin_rho, vmax = vmax_rho)
        colorbar1.update_normal(plot1[0]) # to update the colorbar at each frame
        axis1.set_title(f"Density $p$ at t = {round(frame_number*tstep_diag* dt, 6)} s.")


        vmin_phi = min(min(l) for l in zarray[1,:,2])
        vmax_phi = max(max(l) for l in zarray[1,:,2])

        plot2[0].remove()
        plot2[0] = axis2.plot_trisurf(zarray[1,frame_number,0], zarray[1,frame_number,1], zarray[1,frame_number,2],
                                  cmap = my_cmap, linewidth=0, antialiased=False,
                                  vmin = vmin_phi, vmax = vmax_phi)
        colorbar2.update_normal(plot2[0]) # to update the colorbar at each frame
        axis2.set_title(f"Electrical potential $\\phi$ at t = {round(frame_number*tstep_diag* dt, 6)} s.")

    return update_plot



plot1 = [ax1.plot_trisurf(zarray[0,0,0], zarray[0,0,1], zarray[0,0,2], cmap = my_cmap, linewidth=0, antialiased=False)]
plot2 = [ax2.plot_trisurf(zarray[1,0,0], zarray[1,0,1], zarray[1,0,2], cmap = my_cmap, linewidth=0, antialiased=False)]

cbar1 = fig.colorbar(plot1[0], ax = ax1, shrink = 0.5, aspect = 5)
cbar2 = fig.colorbar(plot2[0], ax = ax2, shrink = 0.5, aspect = 5)



# --- plot and save animation
fps = 10 # frame per sec
ani = FuncAnimation(fig, update_plot_function(fig, ax1, ax2), len(zarray[0]), fargs=(zarray, plot1, cbar1, plot2, cbar2), interval=1000/fps)

plt.show()
ani.save(name_animation + '.mp4',writer='ffmpeg',fps=fps)

