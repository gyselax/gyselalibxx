"""
Save and plot the animation of the density and the electrical potential of 
the diocotron or vortex merger simulation.
"""

import argparse
from pathlib import Path

from RTheta_tool_functions import animate


# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent

perturb_defaults = {'dioctotron':False, 'vortex_merger':True}

parser = argparse.ArgumentParser(description="Plot and save the density and the electrical potential solutions of a given output folder.")

parser.add_argument('-simulation', type=str,
                     choices=['diocotron', 'vortex_merger'])

parser.add_argument('--name', type=str,
                    default = None,
                    help="Name of the saved animation.")

parser.add_argument('--folder', metavar='folder', type=str,
                    default = None,
                    help="Path to the output folder.")

parser.add_argument('--perturb', type=bool, default=None, help='Plot and save only the perturbation (f - f_{eq}).')


args = parser.parse_args()

simulation = args.simulation
name_animation = args.name
folder = args.folder
if_perturb = args.perturb


default_output_dir = current_folder.joinpath(f"../../../build/simulations/geometryRTheta/{simulation}/output/").resolve()


if name_animation is None:
    name_animation = simulation


if args.folder is None:
    folder = default_output_dir
else:
    folder = Path(folder)


if if_perturb is None:
    if_perturb = perturb_defaults[simulation]

# ---------------------------------------------------------------------------


animate(folder, name_animation, if_perturb)
