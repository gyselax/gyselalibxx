"""
Define all the functions needed in the test_convergence.py and 
display_errors.py files. 
"""

import subprocess


# Definition of functions ---------------------------------------
def take_errors_from_out(var_out, degree):
    """
    Take from the output the maximum absolute error 
    for each function.    
    
    Parameters
    ----------
    var_out: string
        Output of the executable tested.
    degree: int 
        Degree of splines influencing the number of
        test functions.

    Returns
    -------
    total_Errors: list of float
        List of errors.
    """
    out_lines = var_out.split('\n')
    out_words = [l.split(' ') for l in out_lines[14+degree:]]
    total_Errors = [float(l[-3]) for l in out_words if len(l)> 2 and l[-3] and l[3] == "Test"]
    return total_Errors


def take_names_from_out(var_out):
    """
    Take from the output the name of each function
    and the degree of splines tested.    
    
    Parameters
    ----------
    var_out : string
        Output of the executable tested.

    Returns
    -------
    degree: int
        Degree of bsplines.
    names: list of string
        List of names.
    """
    out_lines = var_out.split('\n')
    out_words = [l.split(' ') for l in out_lines[:]]
    out_words_flat = [w for l in out_words for w in l]
    degree = int(out_words_flat[out_words_flat.index("degree") + 1])
    names = [l[7] for l in out_words[14+degree:] if len(l)>=7 and l[3] == "Test"]
    return degree, names


def launch_executable(executable, NN):
    """
    Launch executable for grid size of NxN for each N in NN.    
    
    Parameters
    ----------
    executable : string
        Path to the executable.
        
    NN: list of int
        Grid size. 

    Returns
    -------
    total_Errors: list of float
        List of errors computed with take_errors_from_out function.
    names: list of string
        List of names of each test functions not supposed to be exact
        from take_names_from_out function.
    degree: int
        Degree of bsplines from take_names_from_out function.
    """

    N = NN[0]
    with open("grid_size.yaml", "w", encoding="utf-8") as f:
        print("Mesh:", file=f)
        print(f"  r_size: {N}", file=f)
        print(f"  theta_size: {N}", file=f)

    with subprocess.Popen([executable, "grid_size.yaml"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as p:
        out, err = p.communicate()
        assert p.returncode == 0

    print(out)
    degree, names = take_names_from_out(out)
    total_Errors_0 = take_errors_from_out(out, degree)
    if err:
        print(err)
    print("")


    total_Errors = [total_Errors_0]
    for N in NN[1:]:
        with open("grid_size.yaml", "w", encoding="utf-8") as f:
            print("Mesh:", file=f)
            print(f"  r_size: {N}", file=f)
            print(f"  theta_size: {N}", file=f)

        with subprocess.Popen([executable, "grid_size.yaml"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as p:
            out, err = p.communicate()
            assert p.returncode == 0

        print(out)
        total_Errors += [take_errors_from_out(out, degree)]
        if err:
            print(err)
        print("")

    return total_Errors, names, degree
