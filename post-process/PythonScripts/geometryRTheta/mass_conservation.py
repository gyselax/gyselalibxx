"""
Plot the relative variation of the mass of the diocotron or vortex merger solution.
"""

import argparse
from pathlib import Path

from RTheta_tool_functions import plot_mass_conservation


# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent

parser = argparse.ArgumentParser(description="Plot the time evolution of the relative errors on the total mass.")

parser.add_argument('-simulation', type=str,
                     choices=['diocotron', 'vortex_merger'])

parser.add_argument('--name', type=str,
                    default = None,
                    help="Name of the saved animation.")

parser.add_argument('--folder', metavar='folder', type=str,
                    default = None,
                    help="Path to the output folder.")

parser.add_argument('--save', action='store_true', help='Save the plot with the name given by the -name parameter.')


args = parser.parse_args()

simulation = args.simulation
name_file = args.name
folder = args.folder
if_save = args.save


default_output_dir = current_folder.joinpath(f"../../../build/simulations/geometryRTheta/{simulation}/output/").resolve()


if name_file is None:
    name_animation = "mass_conservation_" + simulation


if args.folder is None:
    folder = default_output_dir
else:
    folder = Path(folder)

# ---------------------------------------------------------------------------

plot_mass_conservation(folder, if_save, name_file)
