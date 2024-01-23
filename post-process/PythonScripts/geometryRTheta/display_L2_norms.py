"""
Plot the L2 norms of the perturbed density and the perturbed electrical potential
of the diocotron instabilities simulation.
"""

import argparse
from pathlib import Path

from RTheta_tool_functions import compute_L2_norms

# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent

parser = argparse.ArgumentParser(description="Plot the L2 norms of (phi - phi_{eq}) and (rho - rho_{eq}).")

parser.add_argument('-simulation', type=str,
                     choices=['diocotron'])

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
    name_animation = "L2_norms_" + simulation


if args.folder is None:
    folder = default_output_dir
else:
    folder = Path(folder)
# ---------------------------------------------------------------------------

compute_L2_norms(folder, if_save, name_file)
