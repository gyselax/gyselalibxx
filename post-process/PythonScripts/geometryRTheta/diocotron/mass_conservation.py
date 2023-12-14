"""
Plot the relative variation of the mass of the diocotron solution.
"""

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from gysdata import DiskStore

# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent
default_output_dir = current_folder.joinpath("../../../../build/simulations/geometryRTheta/diocotron/output/").resolve()

parser = argparse.ArgumentParser(description='Plot the time evolution of the relative errors on the total mass.')

parser.add_argument('--folder', metavar='folder', type=str,
                    default = default_output_dir,
                    help=f"Path to the output folder. By default folder = '{default_output_dir}'")

parser.add_argument('-name', metavar='name', type=str, nargs= '?',
                    default = "mass_conservation",
                    help="Name of the saved plot. By default name = 'mass_conservation'.")

parser.add_argument('--save', action='store_true', help='Save the plot with the name given by the -name parameter.')

args = parser.parse_args()

folder = Path(args.folder)
name_file = args.name
if_save = args.save

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

# Treatment for periodicity:
jacobian = np.zeros((Nr+3, Nt+1))
jacobian[:, :-1] = np.array(ds["jacobian"])
jacobian[:, -1] = jacobian[:, 0]



# Compute norms for each time step
Time = np.array(ds["density"].coords["time"])
Mass = np.zeros(iter_nb // tstep_diag + 1)
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

t_step = 10 # step for the time ticks


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

