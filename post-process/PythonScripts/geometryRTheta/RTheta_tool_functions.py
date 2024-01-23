"""
Save and plot the animation of the density and the electrical potential of the vortex merger simulation.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gysdata import DiskStore


# ---------------------------------------------------------------------------

def animate(folder, name_animation, if_perturb):
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
    if_perturb: boolean 
        Indicate if we want to plot the full functions (False) 
        or only the perturbations (True).
    """

    path_data_structure = Path('data_structure_RTheta.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)


    # Get initial data
    dt = float(ds["delta_t"])
    T = float(ds["final_T"])
    iter_nb = int( T / dt ) +1
    tstep_diag = int(ds["time_step_diag"])

    Nr = int(ds["r_size"])
    Nt = int(ds["p_size"])


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
            if if_perturb:
                v_rho = max(-min(min(l) for l in zarray[0,:,2]), max(max(l) for l in zarray[0,:,2]))
                vmin_rho = -v_rho
                vmax_rho = v_rho
            else:
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






def plot_mass_conservation(folder, if_save, name_file):
    """
    Plot the relative variation of the mass of the simulation solution.

    Parameters
    ----------
    folder: Path
        The folder to the output folder containing 
        the output files .h5 of the simulation. 
    if_save: boolean
        If True, save the plot in a file. 
        Otherwise, it doesn't save. 
     name_file: string
        The name of the file where the plot is savec. 
    """

    path_data_structure = Path('data_structure_RTheta.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)

    # Get initial data
    dt = float(ds["delta_t"])
    T = float(ds["final_T"])
    iter_nb = int( T / dt ) +1
    tstep_diag = int(ds["time_step_diag"])

    Nr = int(ds["r_size"])
    Nt = int(ds["p_size"])

    # Treatment for periodicity:
    jacobian = np.zeros((Nr+3, Nt+1))
    jacobian[:, :-1] = np.array(ds["jacobian"])
    jacobian[:, -1] = jacobian[:, 0]



    # Compute norms for each time step
    Time = np.array(ds["density"].coords["time"])
    r_grid = np.array(ds["density"].coords["r"])


    # Treatment for periodicity:
    p_grid = np.zeros(Nt+1)
    p_grid[:-1] = np.array(ds["density"].coords["p"])
    p_grid[-1] = p_grid[0]  + 2*np.pi


    # Get the data at each time step + teatment for periodicity
    rho = np.zeros((iter_nb//tstep_diag+1, Nr+3, Nt+1))
    rho[:, :, :-1] = np.array(ds['density']) [:iter_nb//tstep_diag+1, :, :]
    rho[:, :, -1] = rho[:, :, 0]

    # Get initial mass
    rho_0 = rho[0]
    mass_0 = np.trapz(np.trapz(rho_0 * np.abs(jacobian), p_grid), r_grid)


    absolute_mass = np.trapz(np.trapz(rho * np.abs(jacobian), p_grid), r_grid)
    Mass = (mass_0 - absolute_mass ) / abs(mass_0)



    t_step = 2 # step for the time ticks


    plt.figure(figsize=(10,5))
    plt.title(f"Relative mass variation in time on $N_r \\times N_\\theta =$ {Nr}x{Nt} grid for dt = {dt}.")
    plt.plot(Time, Mass, "-", label="$\\delta \\mathcal{M}(t)$ computed\n $max(\\delta \\mathcal{M}(t)) = $"
            + f"{max(Mass)}"+ "\n $min(\\delta \\mathcal{M}(t)) = $"+ f"{min(Mass)}")

    plt.xlabel(f"$t \\in [{Time[0]}, {Time[-1]}]$ s,  ({round(Time[-1]/dt)} iterations).")
    plt.ylabel("Mass conservation : $\\delta \\mathcal{M}(t) = \\int (p - p_0) / |\\int p_0 |$")
    plt.xticks([t_step*i for i in range(int(T / t_step)+1)])
    plt.legend()
    plt.grid()

    plt.show()

    if if_save:
        plt.savefig(name_file + ".png")







def compute_L2_norms(folder, if_save, name_file):
    """
    Plot the L2 norms of the perturbed density and the perturbed electrical potential
    of the diocotron instabilities simulation.

    Parameters
    ----------
    folder: Path
        The folder to the output folder containing 
        the output files .h5 of the simulation. 
    if_save: boolean
        If True, save the plot in a file. 
        Otherwise, it doesn't save. 
     name_file: string
        The name of the file where the plot is savec. 
    """

    path_data_structure = Path('data_structure_RTheta.yaml')
    ds = DiskStore(folder, data_structure=path_data_structure)


    # Get initial data
    rho_eq = np.array(ds['density_eq'])
    phi_eq = np.array(ds['electrical_potential_eq'])

    jacobian = np.array(ds["jacobian"])

    dt = float(ds["delta_t"])
    T = float(ds["final_T"])
    iter_nb = int( T / dt ) +1
    tstep_diag = int(ds["time_step_diag"])

    Nr = int(ds["r_size"])
    Nt = int(ds["p_size"])
    omega_Im = float(ds["slope"])


    # Compute norms for each time step
    L2norms_rho = np.zeros(iter_nb//tstep_diag+1)
    L2norms_phi = np.zeros(iter_nb//tstep_diag+1)


    # Get the data at each time step
    Time = np.array(ds["density"].coords["time"])
    rho = np.array(ds['density'])
    phi = np.array(ds['electrical_potential'])

    # Compute norms
    L2norms_rho = np.linalg.norm((rho - rho_eq[None,:,:])*abs(jacobian[None,:,:]), 2, axis=(1,2))
    L2norms_phi = np.linalg.norm((phi - phi_eq[None,:,:])*abs(jacobian[None,:,:]), 2, axis=(1,2))





    t_step = 10 # step for the time ticks

    idx_origin = int(30/dt)//tstep_diag
    growth_rate_phi = L2norms_phi[idx_origin]* np.exp((Time-Time[idx_origin])*omega_Im)
    growth_rate_rho = L2norms_rho[idx_origin]* np.exp((Time-Time[idx_origin])*omega_Im)


    # Plot L2 norms
    plt.figure(figsize=(20,7))

    ax = plt.subplot(1,2,1)
    plt.title(f"L2 norm of $\\phi$ in time on $N_r \\times N_\\theta =$ {Nr}x{Nt} grid for dt = {dt}.")
    plt.plot(Time[1:],L2norms_phi[1:], "-", label="computed")
    plt.plot(Time[:int(70/dt)],growth_rate_phi[:int(70/dt)], "--", label = f"growth rate = {omega_Im}")

    ax.set_yscale('log')
    plt.xlabel(f"$t \\in [{Time[0]}, {Time[-1]}]$ s,  ({Time[-1]/dt} iterations).")
    plt.ylabel("$L_2$ norm: $||\\phi-\\phi_{eq}||_{L_2}$")
    plt.xticks([t_step*i for i in range(int(T/t_step)+1)])
    plt.legend()
    plt.grid()


    ax = plt.subplot(1,2,2)
    plt.title(f"L2 norm of $p$ in time on  $N_r \\times N_\\theta =$ {Nr}x{Nt} grid for dt = {dt}")
    plt.plot(Time[1:],L2norms_rho[1:], "-", label="computed")
    plt.plot(Time[:int(70/dt)],growth_rate_rho[:int(70/dt)], "--", label = f"growth rate = {omega_Im}")

    ax.set_yscale('log')
    plt.xlabel(f"$t \\in [{Time[0]}, {Time[-1]}]$ s,  ({round(Time[-1]/dt)} iterations).")
    plt.ylabel("$L_2$ norm: $||p-p_{eq}||_{L_2}$")
    plt.xticks([t_step*i for i in range(int(T/t_step)+1)])
    plt.legend()
    plt.grid()

    plt.show()

    if if_save:
        plt.savefig(name_file + ".png")

