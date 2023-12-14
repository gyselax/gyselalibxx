"""
Plot the L2 norms of the perturbed density and the perturbed electrical potential.
"""

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from gysdata import DiskStore

# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent
default_output_dir = current_folder.joinpath("../../../../build/simulations/geometryRTheta/diocotron/output/").resolve()

parser = argparse.ArgumentParser(description='Plot the L2 norms of (phi - phi_{eq}) and (rho - rho_{eq}).')

parser.add_argument('--folder', metavar='folder', type=str,
                    default = default_output_dir,
                    help=f"Path to the output folder. By default folder = '{default_output_dir}'")

parser.add_argument('-name', metavar='name', type=str, nargs= '?',
                    default = "L2_norms",
                    help="Name of the saved plot. By default name = 'L2_norms'.")

parser.add_argument('--save', action='store_true', help='Save the plot with the name given by the -name parameter.')

args = parser.parse_args()

folder = Path(args.folder)
name_file = args.name
if_save = args.save

# ---------------------------------------------------------------------------

path_data_structure = Path('data_structure_diocotron.yaml')
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

