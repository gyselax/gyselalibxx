"""
Define all the functions needed in the python files of the advection_2d_rp folder.
"""

import subprocess

import argparse

import numpy as np



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

    print("")

    args = parser.parse_args()

    rmin = args.rmin
    rmax = args.rmax
    Nr = args.Nr
    Nth = args.Nth
    dt = args.dt
    T = args.T
    curves = args.curves
    feet = args.feet
    if_plot = args.plot

    executable = args.executable[0]

    return executable, rmin, rmax, Nr, Nth, dt, T, curves, feet, if_plot



def params_file(rmin, rmax, Nr, Nth, dt, T, curves=False, feet=False):
    """
    Fill in the params.yaml parameters file with the correct inputs.

    Parameters
    ----------
    rmin : float
        Minimum value of the r values.

    rmax : float
        Maximum value of the r values.

    Nr : int
        Number of break points in r dimension

    Nth : int
        Number of break points in theta dimension

    dt : float
        Time step.

    T : float
        Final time.

    curves : bool
        Boolean to select if the values of the advected function are saved.

    feet : bool
        Boolean to select if the characteristic feet are saved.
    """
    with open("params.yaml", "w", encoding="utf-8") as f:
        print("Mesh:", file=f)
        print(f"  r_min: {rmin}", file=f)
        print(f"  r_max: {rmax}", file=f)
        print(f"  r_size: {Nr}", file=f)
        print(f"  p_size: {Nth}", file=f)
        print(f"  time_step: {dt}", file=f)
        print(f"  final_time: {T}", file=f)
        print("  save_curves: "+"true"*curves + "false"*(not curves), file=f)
        print("  save_feet: "+"true"*feet + "false"*(not feet), file=f)


def execute(executable, rmin, rmax, Nr, Nth, dt, T, curves=False, feet=False, print_out=True):
    """
    Launch the executable given as input.

    Parameters
    ----------
    executable : string
        Path the executable of the test we want to launch.

    rmin : float
        Minimum value of the r values.

    rmax : float
        Maximum value of the r values.

    Nr : int
        Number of break points in r dimension

    Nth : int
        Number of break points in theta dimension

    dt : float
        Time step.

    T : float
        Final time.

    curves : bool
        Boolean to select if the values of the advected function are saved.

    feet : bool
        Boolean to select if the characteristic feet are saved.

    print_out : bool
        Boolean to select if the output are printed in the console.

    Returns
    -------
        The output of the executable.
    """
    params_file(rmin, rmax, Nr, Nth, dt, T, curves, feet)

    with subprocess.Popen([executable, "params.yaml"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as p:
        out, err = p.communicate()
        assert p.returncode == 0

    if print_out:
        print(out)
    if err:
        print(err)
    if print_out:
        print("")

    return out



def get_fct_names(var_out):
    """
    Get the reduced names of all the test cases from the output.
    The names are reduced to the name of the simulation (no mapping,
    no advection domain, no method).

    Parameters
    ----------
    var_out : string
        The output of the executable.

    Returns
    -------
        A list of string containing the names of the test cases.
    """
    out_lines = var_out.split('\n')
    out_words = [line.split(' - ') for line in out_lines[1:] if "MAPPING" in line and "DOMAIN" in line]
    fct_names = [line[3][:-3].lower() for line in out_words]
    return fct_names

def get_full_fct_names(var_out):
    """
    Get the full names of all the test cases from the output.

    Parameters
    ----------
    var_out : string
        The output of the executable.

    Returns
    -------
        A list of string containing the names of the test cases.
    """
    out_lines = var_out.split('\n')
    out_words = [[l.strip(' :') for l in line.split(' - ')] for line in out_lines[1:] if "MAPPING" in line and "DOMAIN" in line]
    fct_names = [line[3].capitalize() + " with " + line[2].lower() + " on "
                 + line[0].lower() + " and " + line[1].lower()  for line in out_words]
    return fct_names

def get_fct_name_key(full_name):
    """
    Get the keys which identify the test case from its full name.

    Parameters
    ----------
    full_name : str
        The full name of the case as outputted by get_full_fct_names.

    Returns
    -------
    problem_type : str
        The key describing the problem type ['Translation'|'Rotation'|'Decentered rotation'].
    time_integration_method : str
        The time integration method used to solve the problem ['euler', 'crank nicolson', 'rk3', 'rk4'].
    mapping : str
        The mapping which was examined in the test ['circular', 'czarny_physical', 'czarny_pseudo_cart', 'discrete'].
    """
    problem_type, s = full_name.split(' with ')
    time_integration_method, s = s.split(' on ')
    mapping, domain = s.split(' and ')
    mapping, _ = mapping.split(' mapping')
    if mapping == 'czarny':
        domain, _ = domain.split(' domain')
        domain = domain.replace(' ','_').lower()
        mapping = f'{mapping}_{domain}'
    return problem_type, time_integration_method, mapping


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
    File = np.loadtxt(namefile)
    #Ir, Ip, IVr, IVp = [], [], [], []
    CoordR = File[:, 4]
    CoordP = File[:, 5]
    CoordX = File[:, 6]
    CoordY = File[:, 7]
    list_F = File[:, 8]
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
        The distance between the two poins.
    """
    Distance = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    return Distance
