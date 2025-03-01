"""
Plot the L2 norm of the density and the electrostatic potential of
the guiding centre XY simulation.
"""

import argparse
from pathlib import Path

from XY_tool_functions import plot_L2_norm

# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent

parser = argparse.ArgumentParser(description="Plot the L2 norms of the density and the electrostatic potential solutions of a given output folder at all time steps.")

parser.add_argument('--name', type=str,
                    default = None,
                    help="Name of the saved animation. If any name is given, it doesn't save the plot.")

parser.add_argument('--folder', metavar='folder', type=str,
                    default = None,
                    help="Path to the output folder. Default path from current folder: ../../../build/simulations/geometryXY/guiding_center/output/")

args = parser.parse_args()

name_file = args.name
folder = args.folder

default_output_dir = current_folder.joinpath("../../../build/simulations/geometryXY/guiding_center/output/").resolve()

if args.folder is None:
    folder = default_output_dir
else:
    folder = Path(folder)

# ---------------------------------------------------------------------------
plot_L2_norm(folder, name_file)
