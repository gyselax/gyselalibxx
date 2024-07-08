"""
Save and plot the mass conservation of the guiding center XY simulation.
"""

import argparse
from pathlib import Path

from XY_tool_functions import plot_mass_conservation

# ---------------------------------------------------------------------------
current_folder = Path(__file__).parent

parser = argparse.ArgumentParser(description="Plot and save the mass conservation of a given output folder.")

parser.add_argument('--name', type=str,
                    default = None,
                    help="Name of the saved figure. If any name is given, it doesn't save the plot.")

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
plot_mass_conservation(folder, name_file)
