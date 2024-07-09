"""
Define useful functions for the current folder.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gysdata import DiskStore

# FONCTIONS USED IN THIS FILE ####################################################################
def get_data(folder):
    """
    Get data from the output folder.

    Parameters
    ----------
    folder: Path
        The folder to the output folder containing
        the output files .h5 of the simulation.

    Returns
    -------
        A dictionnary containing
            Nx: Number of cells in x.
            Ny: Number of cells in y.
            x_grid: Grid in x.
            x_grid: Grid in y.
            dt: Time step.
            T: Final time.
            tstep_diag: Time step applied to save data.
            Time: Array containing all the saved time steps.
    """
    path_data_structure = Path('data_structure_XY.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get data
    return  {"Nx": int(ds["Nx_spline_cells"]),
            "Ny": int(ds["Ny_spline_cells"]), 
            "x_grid": np.array(ds['fdistribu'].coords['x'].values), 
            "y_grid": np.array(ds['fdistribu'].coords['y'].values), 
            "dt": float(ds["delta_t"]), 
            "T": float(ds["final_time"]), 
            "tstep_diag": int(ds["nbstep_diag"]), 
            "Time": np.array(ds['fdistribu'].coords['time'].values)}


def get_periodic_grid(grid):
    """
    To use the trapezoid method, we need to duplicate
    points to make the grids periodic for the
    np.trapz() operator.
    We supose the grid uniform.

    Parameters
    ----------
    grid: array
        Grid containing the coordinates without duplication.

    Returns
    -------
        Grid with duplication point for applying np.trapz().
    """
    N = len(grid)
    periodic_grid = np.empty(N+1)
    periodic_grid[:-1] = grid
    periodic_grid[-1] = periodic_grid[-2] + (periodic_grid[-2] -periodic_grid[-3])
    return periodic_grid




def get_periodic_function(data):
    """
    To use the trapezoid method, we need to duplicate
    points to make the functions periodic for
    the np.trapz() operator.
    For function defined on space domain or time and space domain.
    The periodicity treatment is applied to the two last 
    dimensions. 

    Parameters
    ----------
    data: array
        Array containing the data without duplication.

    Returns
    -------
        Array with duplication point for applying np.trapz().
    """
    shape = list(data.shape)
    shape[-1] += 1 # Add 1 to Ny
    shape[-2] += 1 # Add 1 to Nx
    periodic_data = np.empty(shape)
    periodic_data[..., :-1, :-1] = data
    periodic_data[..., -1, :-1]  = periodic_data[..., 0, :-1]
    periodic_data[..., :-1, -1]  = periodic_data[..., :-1, 0]
    periodic_data[..., -1,-1] = periodic_data[...,0,0]
    return periodic_data




# FONCTIONS FOR OTHER FILES ######################################################################
def animate(folder, name_animation):
    """
    Save and plot the animation of the density and the electrical potential of the simulation.

    Parameters
    ----------
    folder: Path
        The folder to the output folder containing
        the output files .h5 of the simulation.
    name_animation: string
        The name of the file where the animation is
        saved.
    """
    path_data_structure = Path('data_structure_XY.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get data
    data = get_data(folder)
    Time = data["Time"]
    X, Y = np.meshgrid(data["x_grid"], data["y_grid"])

    # --- function data
    F = np.array(ds['fdistribu'])
    Phi = np.array(ds['electrostatic_potential'])

    zarray = np.array([F, Phi])

    # --- min and max values for colorbar
    vmin_F = F.min()
    vmax_F = F.max()

    vmin_phi = Phi.min()
    vmax_phi = Phi.max()

    # Animation
    plt.rc('text', usetex=True) # for latex font on plots.
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16,6))

    for ax in [ax1, ax2]:
        ax.set_aspect("equal")
        ax.set_xlabel(r"x")
        ax.set_ylabel(r"y")

    my_cmap = plt.get_cmap('jet')#('Spectral')#('plasma')#('jet')#('gnuplot')#('seismic')
    def update_plot_function(axis1, axis2):
        """
        Define a update_plot function which update the plots
        in the FuncAnimation function.

        Parameters
        ----------
        axis1: matplotlib.pyplot.axis
            The first subplot object.
        axis2: matplotlib.pyplot.axis
            The second subplot object.

        Returns
        -------
            update_plot, the function which updates the plots
            for the FuncAnimation function.
        """
        def update_plot(frame_number, zarray, plot1, colorbar1, plot2, colorbar2):
            plot1[0].remove()
            plot1[0] = axis1.pcolormesh(X, Y, zarray[0,frame_number].T,
                                      cmap = my_cmap, vmin = vmin_F, vmax = vmax_F)
            colorbar1.update_normal(plot1[0]) # to update the colorbar at each frame
            axis1.set_title(f"Density $f$ at t = {round(Time[frame_number], 6)} s.")


            plot2[0].remove()
            plot2[0] = axis2.pcolormesh(X, Y, zarray[1,frame_number].T,
                                      cmap = my_cmap, vmin = vmin_phi, vmax = vmax_phi)
            colorbar2.update_normal(plot2[0]) # to update the colorbar at each frame
            axis2.set_title(f"Electrical potential $\\phi$ at t = {round(Time[frame_number], 6)} s.")
        return update_plot

    plot1 = [ax1.pcolormesh(X, Y, zarray[0,0].T, cmap = my_cmap, vmin = vmin_F, vmax = vmax_F)]
    plot2 = [ax2.pcolormesh(X, Y, zarray[1,0].T, cmap = my_cmap, vmin = vmin_phi, vmax = vmax_phi)]

    cbar1 = fig.colorbar(plot1[0], ax = ax1, shrink = 0.5, aspect = 5)
    cbar2 = fig.colorbar(plot2[0], ax = ax2, shrink = 0.5, aspect = 5)

    # --- plot and save animation
    fps = 10 # frame per sec
    ani = FuncAnimation(fig, update_plot_function(ax1, ax2), len(zarray[0]), fargs=(zarray, plot1, cbar1, plot2, cbar2), interval=1000/fps)

    plt.show()
    if name_animation is not None:
        ani.save(name_animation + '.mp4',writer='ffmpeg',fps=fps)






def plot_mass_conservation(folder, name_file):
    """
    Plot the variation of the mass of the simulation solution.

    Parameters
    ----------
    folder: Path
        The path to the `output/` folder containing
        the output files .h5 of the simulation.
     name_file: string
        The name of the file where the plot is saved.
    """
    path_data_structure = Path('data_structure_XY.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get initial data
    data = get_data(folder)
    Nx, Ny = data["Nx"], data["Ny"]
    dt, T, Time = data["dt"], data["T"], data["Time"]

    # Treatment for periodicity
    x_grid = get_periodic_grid(data["x_grid"])
    y_grid = get_periodic_grid(data["y_grid"])

    # --- get initial mass
    fdistribu_init = get_periodic_function(ds['fdistribu'][0])
    mass_init = np.trapz(np.trapz(fdistribu_init, y_grid, axis=1), x_grid)

    # --- get mass at each time step
    fdistribu = get_periodic_function(ds['fdistribu'])
    absolute_mass = np.trapz(np.trapz(fdistribu, y_grid, axis=2), x_grid, axis=1)

    # --- get absolutse mass evolution
    Mass = absolute_mass - mass_init


    t_step =  T//10 # step for the time ticks

    plt.figure(figsize=(10,5))
    plt.rc('text', usetex=True)
    plt.title(f"Mass variation in time on $N_x \\times N_y =$ {Nx}x{Ny} grid for dt = {dt}.")
    plt.plot(Time, Mass, "-", label="$\\delta \\mathcal{M}(t)$ computed\n $max(\\delta \\mathcal{M}(t)) = $"
            + f"{max(Mass):.6e}"+ "\n $min(\\delta \\mathcal{M}(t)) = $"+ f"{min(Mass):.6e}")

    plt.xlabel(f"$t \\in [{Time[0]}, {Time[-1]}]$ s,  ({round(Time[-1]/dt)} iterations).")
    plt.ylabel(r"Mass conservation : $\delta \mathcal{M}(t) = \int f - \int f_0$")
    plt.xticks([t_step*i for i in range(int(T / t_step)+1)])
    plt.legend()
    plt.grid()
    plt.show()

    if name_file is not None:
        plt.savefig(name_file + ".png")






def plot_energy_conservation(folder, name_file):
    """
    Plot the variation of the energie of the simulation solution.

    Parameters
    ----------
    folder: Path
        The path to the `output/` folder containing
        the output files .h5 of the simulation.
     name_file: string
        The name of the file where the plot is saved.
    """
    path_data_structure = Path('data_structure_XY.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get initial data
    data = get_data(folder)
    Nx, Ny = data["Nx"], data["Ny"]
    dt, T, Time = data["dt"], data["T"], data["Time"]


    # Treatment for periodicity
    x_grid = get_periodic_grid(data["x_grid"])
    y_grid = get_periodic_grid(data["y_grid"])


    # Get data
    # --- get initial energy
    electric_field_x_init = ds['electric_field_x'][0]
    electric_field_y_init = ds['electric_field_y'][0]

    norm_electric_field_init = get_periodic_function(electric_field_x_init**2 + electric_field_y_init**2)
    energy_init = np.trapz(np.trapz(norm_electric_field_init, y_grid, axis =1), x_grid)

    # --- get energy for each time step
    electric_field_x = ds['electric_field_x']
    electric_field_y = ds['electric_field_y']

    norm_electric_field = get_periodic_function(electric_field_x**2 + electric_field_y**2)
    absolute_energy = np.trapz(np.trapz(norm_electric_field, y_grid, axis = 2), x_grid, axis = 1)

    # --- variation energy
    Energy = absolute_energy - energy_init


    t_step =  T//10 # step for the time ticks

    plt.figure(figsize=(10,5))
    plt.rc('text', usetex=True)
    plt.title(f"Energy variation in time on $N_x \\times N_y =$ {Nx}x{Ny} grid for dt = {dt}.")
    plt.plot(Time, Energy, "-", label="$\\delta \\mathcal{W}(t)$ computed\n $max(\\delta \\mathcal{W}(t)) = $"
            + f"{max(Energy):.6e}"+ "\n $min(\\delta \\mathcal{W}(t)) = $"+ f"{min(Energy):.6e}")

    plt.xlabel(f"$t \\in [{Time[0]}, {Time[-1]}]$ s,  ({round(Time[-1]/dt)} iterations).")
    plt.ylabel(r"Energie conservation : $\delta \mathcal{W}(t) = \int |E|^2 - \int |E_0|^2$")
    plt.xticks([t_step*i for i in range(int(T / t_step)+1)])
    plt.legend()
    plt.grid()
    plt.show()

    if name_file is not None:
        plt.savefig(name_file + ".png")






def plot_L2_norm(folder, name_file):
    """
    Plot the L2 norm of the density and the electrostatic potential.

    Parameters
    ----------
    folder: Path
        The path to the `output/` folder containing
        the output files .h5 of the simulation.
     name_file: string
        The name of the file where the plot is saved.
    """
    path_data_structure = Path('data_structure_XY.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get initial data
    data = get_data(folder)
    Nx, Ny = data["Nx"], data["Ny"]
    dt, T, Time = data["dt"], data["T"], data["Time"]

    # --- get equilibrium data
    fdistribu_eq = np.array(ds['fdistribu_equilibrium'])

    # --- get data at each time step
    fdistribu = np.array(ds['fdistribu'])
    electrostatic_potential = np.array(ds['electrostatic_potential'])

    L2norm_F = np.linalg.norm((fdistribu - fdistribu_eq[None, :, :]), 2, axis=(1,2))
    L2norm_Phi = np.linalg.norm((electrostatic_potential), 2, axis=(1,2))


    t_step = T//10 # step for the time ticks

    plt.figure(figsize=(20,7))
    plt.rc('text', usetex=True)

    ax = plt.subplot(1,2,1)
    plt.title(r"L2 norm $||f - f_{eq}||_{\mathcal{L}^2}$ in time on " + f"$N_x \\times N_y =$ {Nx}x{Ny} grid for dt = {dt}.")
    plt.plot(Time, L2norm_F, "-", label="computed")
    ax.set_yscale('log')
    plt.xlabel(f"$t \\in [{Time[0]}, {Time[-1]}]$ s,  ({Time[-1]/dt} iterations).")
    plt.ylabel(r"L2 norm $||f -f_{eq}||_{\mathcal{L}^2}$")
    plt.xticks([t_step*i for i in range(int(T/t_step)+1)])
    plt.legend()
    plt.grid()

    ax = plt.subplot(1,2,2)
    plt.title(r"L2 norm $||\phi||_{\mathcal{L}^2}$ in time on " + f"$N_x \\times N_y =$ {Nx}x{Ny} grid for dt = {dt}.")
    plt.plot(Time, L2norm_Phi, "-", label="computed")
    ax.set_yscale('log')
    plt.xlabel(f"$t \\in [{Time[0]}, {Time[-1]}]$ s,  ({Time[-1]/dt} iterations).")
    plt.ylabel(r"L2 norm $||\phi||_{\mathcal{L}^2}$")
    plt.xticks([t_step*i for i in range(int(T/t_step)+1)])
    plt.legend()
    plt.grid()

    plt.show()


    if name_file is not None:
        plt.savefig(name_file + ".png")






def plot_individual_curve(folder, time, name_file):
    """
    Save and plot the animation of the density and the electrical potential of the simulation.

    Parameters
    ----------
    folder: Path
        The folder to the output folder containing
        the output files .h5 of the simulation.
    name_animation: string
        The name of the file where the animation is
        saved.
    """
    path_data_structure = Path('data_structure_XY.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get initial data
    data = get_data(folder)
    dt, tstep_diag, Time = data["dt"], data["tstep_diag"], data["Time"]
    X, Y = np.meshgrid(data["x_grid"], data["y_grid"])

    # --- function data
    time_step = int(time / dt)
    idx =  int(time_step / tstep_diag)
    F = np.array(ds['fdistribu'])[idx]
    Phi = np.array(ds['electrostatic_potential'])[idx]


    # Plot curves
    plt.rc('text', usetex=True)
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16,4))
    my_cmap = plt.get_cmap('jet')#('Spectral')#('plasma')#('jet')#('gnuplot')

    vmin_F = F.min()
    vmax_F = F.max()
    plot1 = ax1.pcolormesh(X, Y, F.T, cmap = my_cmap, vmax= vmax_F, vmin=vmin_F)
    fig.colorbar(plot1, ax = ax1)
    ax1.set_title(r"Density $f$ "+f"at t = {round(Time[idx], 6)} s.")

    vmin_phi = Phi.min()
    vmax_phi = Phi.max()
    plot2 = ax2.pcolormesh(X, Y, Phi.T, cmap = my_cmap, vmax= vmax_phi, vmin=vmin_phi)
    fig.colorbar(plot2, ax = ax2)
    ax2.set_title(r"Electrical potential $\phi$ "+f"at t = {round(Time[idx], 6)} s.")

    for ax in [ax1, ax2]:
        ax.set_aspect("equal")
        ax.set_xlabel(r"x")
        ax.set_ylabel(r"y")
    plt.show()

    if name_file is not None:
        plt.savefig(name_file + ".png")
