"""
Define all the functions needed in the python files of the advection_2d_rp folder.
"""

import argparse
import os
import subprocess

from collections import namedtuple

import numpy as np
import yaml

SimulationResults = namedtuple('SimulationResults', ('name', 'l_inf_error', 'l_2_error'))

def set_input(rmin_def, rmax_def, Nr_def, Nth_def, dt_def, T_def, curves_def=False, feet_def=False):
    """
    Use the argparse package to fill in the params.yaml input parameters file.
    Possibility to change the values of the input parameters in the commande line.

    Parameters
    ----------
    rmin_def : float
        Minimum value of the r values.

    rmax_def : float
        Maximum value of the r values.

    Nr_def : int
        Number of break points in r dimension

    Nth_def : int
        Number of break points in theta dimension

    dt_def : float
        Time step.

    T_def : float
        Final time.

    curves_def : bool
        Boolean to select if the values of the advected function are saved.

    feet_def : bool
        Boolean to select if the characteristic feet are saved.


    Returns
    -------
        The executable, rmin, rmax, Nr, Nth, dt, T, curves, feet parameters
        for the simulation. It can by the default ones or the modified ones.
    """
    # Inputs --------------------------------------------------------
    parser = argparse.ArgumentParser(description='Launch the executable given with the following parameters.')
    parser.add_argument('executable', metavar='executable', type=str, nargs= 1, help='Path to be executable.')
    parser.add_argument('-rmin', metavar='rmin', type=float, nargs='?', default = rmin_def,
                        help=f'Minimum value for r. By default, rmin = {rmin_def}.')
    parser.add_argument('-rmax', metavar='rmax', type=float, nargs='?', default = rmax_def,
                        help=f'Maximum value for r. By default, rmax = {rmax_def}.')

    parser.add_argument('-Nr', metavar='Nr', type=int, nargs='?', default = Nr_def,
                        help=f'Number of break points in r dimension. By default, Nr = {Nr_def}.')
    parser.add_argument('-Nth', metavar='Nth', type=int, nargs='?', default = Nth_def,
                        help=f'Number of break points in theta dimension. By default, Nth = {Nth_def}.')

    parser.add_argument('-dt', metavar='dt', type=float, nargs='?', default = dt_def,
                        help=f'Time step. By default, dt = {dt_def}.')
    parser.add_argument('-T', metavar='T', type=float, nargs='?', default = T_def,
                        help=f'Final time. By default, dt = {T_def}.')

    parser.add_argument('-curves', metavar='curves', type=bool, nargs='?', default = curves_def,
                        help=f'Boolean to select if the values of the advected function are saved. By default, curves = {"True"*curves_def + "False"*(not curves_def)}.')
    parser.add_argument('-feet', metavar='feet', type=bool, nargs='?', default = feet_def,
                        help=f'Boolean to select if the feet are saved. By default, feet = {"True"*feet_def + "False"*(not feet_def)}.')

    parser.add_argument('--plot', action='store_true', help='Plot the results.')
    parser.add_argument('--verbose', action='store_true', help='Output information about the execution.')

    args = parser.parse_args()

    yaml_parameters = {'SplineMesh': {
                            'r_min': args.rmin,
                            'r_max': args.rmax,
                            'r_ncells': args.Nr,
                            'p_ncells': args.Nth
                        },
                        'Time': {
                            'time_step': args.dt,
                            'final_time': args.T
                        },
                        'Output': {
                            'save_curves': args.curves,
                            'save_feet': args.feet
                        }
            }

    executable = args.executable[0]

    return executable, yaml_parameters, args.plot, args.verbose


def execute(executable, yaml_parameters, print_out=True):
    """
    Launch the executable given as input.

    Parameters
    ----------
    executable : string
        Path the executable of the test we want to launch.

    yaml_parameters : dict
        A dictionary describing the simulation parameters.

    print_out : bool
        Boolean to select if the output are printed in the console.

    Returns
    -------
        The output of the executable.
    """
    params_file = executable+"_params.yaml"
    with open(params_file, "w", encoding="utf-8") as f:
        print(yaml.dump(yaml_parameters), file=f)

    with subprocess.Popen([executable, params_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as p:
        out, err = p.communicate()
        if err:
            print(err)
        assert p.returncode == 0

    if print_out:
        print(out)
        print("")

    return out

def get_simulation_config(executable_file):
    """
    Get the keys which identify the simulation from the executable name.

    Parameters
    ----------
    executable_file : str
        The location of the executable.

    Returns
    -------
    mapping : str
        The key which identifies the mapping (and domain).
    method : str
        The key which identifies the numerical method.
    domain : str
        The key which identifies the domain (physical/pseudo cartesian).
    simulation : str
        The key which identifies the simulation.
    description : str
        The description of the simulation.
    """
    executable_name = os.path.basename(executable_file)
    mapping, method, simulation = executable_name.replace('advection_', '').split('__')
    mapping, domain = mapping.lower().split('_mapping_')
    method = method.replace('_METHOD', '').lower()
    simulation = simulation.replace('_SIMULATION', '').capitalise()
    description = f'{simulation.capitalise()} with {method} on {mapping} and {domain}'
    return mapping, method, domain, simulation, description

def extract_simulation_results(var_out):
    """
    Get the results of a set of simulations.

    Parameters
    ----------
    var_out : list[str]
        A list containing the output of each simulation that was run.

    Returns
    -------
    dict[tuple[str, str, str], SimulationResults]
        A dictionary whose keys are a 3-tuple describing the numerical method, the mapping, and the
        simulation that was run, and whose values are SimulationResults tuples containing the values
        for the different errors and a descriptive name for the simulation.
    """
    out_lines = [o.split('\n') for o in var_out]

    # Extract lines which describe the simulation
    # E.g:     CIRCULAR MAPPING - PHYSICAL DOMAIN - CRANK NICOLSON - DECENTRED ROTATION :
    simulation_description = [[l.strip(' :').lower() for l in line.split(' - ')] for line in out_lines[0]
                                         if "MAPPING" in line and "DOMAIN" in line]

    # Build a readable description of the simulation
    simulation_names = [f'{simu.capitalise()} with {method} on {mapping} and {domain}'
                        for mapping, domain, method, simu in simulation_description]

    simulation_keys = [(method, mapping.replace(" mapping","")+" "+domain.replace(" domain",""), simu.capitalise())
                        for mapping, domain, method, simu in simulation_description]

    l_inf_errors = [[float(line.split()[-1]) for line in o if "Max absolute error" in line] for o in out_lines]
    l_2_errors = [[float(line.split()[-1]) for line in o if "Relative L2 norm error" in line] for o in out_lines]

    simulation_results = {key: SimulationResults(name, linf_err, l2_err)
                  for key, name, linf_err, l2_err
                    in zip(simulation_keys, simulation_names, zip(*l_inf_errors), zip(*l_2_errors))}

    return simulation_results


def treatment(namefile):
    """
    Get lists of function data from a text file.

    Parameters
    ----------
    namefile : string
        Path to the data file.

    Returns
    -------
        Lists of the function values; coordinates in r;
        coordinates in theta; coordinates in x; coordinates in y.
    """
    File = np.loadtxt(namefile, float)
    #Ir, Ip, IVr, IVtheta = [], [], [], []
    CoordR = File[:, 2]
    CoordP = File[:, 3]
    CoordX = File[:, 4]
    CoordY = File[:, 5]
    list_F = File[:, 6]
    return list_F, CoordR, CoordP,  CoordX, CoordY




def treatment_feet(namefile):
    """
    Get lists of characteristic feet data from a text file.

    Parameters
    ----------
    namefile : string
        Path to the data file.

    Returns
    -------
        Lists of the index in r; index in theta; coordinates in r;
        coordinates in theta; coordinates in x; coordinates in y;
        feet in r; feet in theta; feet in x; feet in y.
    """
    Ir, Ip, Cr, Cp, Cx, Cy, Fr, Fp, Fx, Fy = np.loadtxt(namefile).T
    return np.int64(Ir), np.int64(Ip), Cr, Cp, Cx, Cy, Fr, Fp, Fx, Fy


def take_errors_from_out(var_out):
    """
    Get the errors from the output of the executable.

    Parameters
    ----------
    var_out : string
        The output of the executable.

    Returns
    -------
        List containing the maximum absolute errors given in the
        output of the executable for each functions.
    """
    out_lines = var_out.split('\n')
    out_words = [line.split(' ') for line in out_lines[1:]]
    total_Errors = [float(line[-1]) for line in out_words if len(line)>2 and line[1] == "-"]
    return total_Errors


def distance(x1, y1, x2, y2):
    """
    Compute the distance between two points in the
    cartesian domain.

    Parameters
    ----------
    x1 : array
        First coordinates of the first point.
    y1 : array
        Second coordinates of the first point.

    x2 : array
        First coordinates of the second point.
    y2 : array
        Second coordinates of the second point.

    Returns
    -------
        The distance between the two points.
    """
    Distance = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    return Distance
