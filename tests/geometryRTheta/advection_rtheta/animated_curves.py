"""
Save an animated version of the selected advection simulations.

Parameters
----------
executable : string
    Path to the executable of the advection simulation.

rmin : float
    Minimum value of the r values.
    Possibility for the user to change it by adding the two
    new sizes in the command line.

rmax : float
    Maximum value of the r values.
    Possibility for the user to change it by adding the two
    new sizes in the command line.

Nr : int
    Number of break points in r dimension
    Possibility for the user to change it by adding the two
    new sizes in the command line.

Nth : int
    Number of break points in theta dimension
    Possibility for the user to change it by adding the two
    new sizes in the command line.

dt : float
    Time step.
    Possibility for the user to change it by adding the two
    new sizes in the command line.

T : float
    Final time.
    Possibility for the user to change it by adding the two
    new sizes in the command line.

curves : bool
    Boolean to select if the values of the advected function are saved.
    Possibility for the user to change it by adding the two
    new sizes in the command line.

feet : bool
    Boolean to select if the characteristic feet are saved.
    Possibility for the user to change it by adding the two
    new sizes in the command line.

Returns
-------
    Save an animation of the selected advection simulations.
"""

import os
import shutil

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation

from advection_functions import set_input, execute, treatment, get_simulation_config



# ---------------------------------------------------------------------------
def animate (iter_nb, file_name, test_name, output_folder):
    """
    Save an animated version of the selected advection simulations.

    Parameters
    ----------
    iter_nb : int
        The number of total iterations in the simulation.
        It is also the number of files with the saved data.

    file_name : string
        The name of the saved file.

    test_name : string
        The title of the test simulation.
    """
    iter_nb_lim = 100
    if iter_nb < iter_nb_lim:
        iter_numbers = range(iter_nb)
    else:
        divisor = iter_nb//iter_nb_lim + 1
        iter_numbers = range(0,iter_nb, divisor)

    zarray = []
    for idx in iter_numbers:
        namefile = os.path.join(output_folder, f"after_{idx}.txt")
        list_F, _, _,  CoordX, CoordY = treatment(namefile)

        zarray += [[CoordX, CoordY, list_F]]

    zarray = np.array(zarray)

    x_max, x_min = max(CoordX), min(CoordX)
    y_max, y_min = max(CoordY), min(CoordY)

    fps = 10 # frame per sec
    my_cmap = plt.get_cmap('inferno')

    def update_plot_function(figure, axis1, axis2):
        def update_plot(frame_number, zarray, plot1, colorbar1, plot2, colorbar2):
            vmin = min(min(l) for l in zarray[:,2])
            vmax = max(max(l) for l in zarray[:,2])

            plot1[0].remove()
            plot1[0] = axis1.plot_trisurf(zarray[frame_number,0], zarray[frame_number,1], zarray[frame_number,2],
                                      cmap = my_cmap, linewidth=0, antialiased=False,
                                      vmin = vmin, vmax = vmax)
            colorbar1.update_normal(plot1[0]) # to update the colorbar at each frame


            plot2[0].remove()
            plot2[0] = axis2.plot_trisurf(zarray[frame_number,0], zarray[frame_number,1], zarray[frame_number,2],
                                      cmap = my_cmap, linewidth=0, antialiased=False,
                                      vmin = vmin, vmax = vmax)
            colorbar2.update_normal(plot2[0]) # to update the colorbar at each frame

        return update_plot

    folder = "animations/"
    if not os.path.exists(folder):
        os.mkdir(folder)


    def set_axis (figure, idx, elev, azim, roll):
        ax = figure.add_subplot(idx, projection='3d')
        ax.view_init(elev, azim, roll)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_xlim(min(x_min,y_min), max(x_max,y_max))
        ax.set_ylim(min(x_min,y_min), max(x_max,y_max))
        ax.set_title(test_name)

        plot = [ax.plot_trisurf(zarray[0,0], zarray[0,1], zarray[0,2], cmap = my_cmap, linewidth=0, antialiased=False)]
        cbar = figure.colorbar(plot[0], ax = ax, shrink = 0.5, aspect = 10)

        return ax, plot, cbar


    fig = plt.figure(figsize=(16,8))

    ax1, plot1, cbar1 = set_axis(fig, 121, 90, -90, 0)
    ax2, plot2, cbar2 = set_axis(fig, 122, 60, -45, 0)

    ani = FuncAnimation(fig, update_plot_function(fig, ax1, ax2), len(zarray), fargs=(zarray, plot1, cbar1, plot2, cbar2), interval=1000/fps)

    #ani.save(folder + file_name+".mp4",writer='ffmpeg',fps=fps)
    ani.save(folder + file_name +".gif",writer='imagemagick',fps=fps)

    plt.clf()


explanations  = """
      Euler                     Crank Nicolson                    RK3                          RK4
 0. translation on A         12. translation on A         24. translation on A         36. translation on A
 1. rotation on A            13. rotation on A            25. rotation on A            37. rotation on A
 2. decentred rotation on A  14. decentred rotation on A  26. decentred rotation on A  38. decentred rotation on A

 3. translation on B         15. translation on B         27. translation on B         39. translation on B
 4. rotation on B            16. rotation on B            28. rotation on B            40. rotation on B
 5. decentred rotation on B  17. decentred rotation on B  29. decentred rotation on B  41. decentred rotation on B

 6. translation on C         18. translation on C         30. translation on C         42. translation on C
 7. rotation on C            19. rotation on C            31. rotation on C            43. rotation on C
 8. decentred rotation on C  20. decentred rotation on C  32. decentred rotation on C  44. decentred rotation on C

 9. translation on D         21. translation on D         33. translation on D         45. translation on D
10. rotation on D            22. rotation on D            34. rotation on D            46. rotation on D
11. decentred rotation on D  23. decentred rotation on D  35. decentred rotation on D  47. decentred rotation on D

 with A : Circular mapping and advection on the physical domain,
      B : Czarny mapping and advection on the physical domain,
      C : Czarny mapping and advection on the pseudo Cartesian domain,
      D : Discrete mapping and advection on the pseudo Cartesian domain.

"""


# Get the inputs -----------------------------------------------
executable, yaml_parameters, _, verbose = set_input(0, 1, 60, 120, 0.01, 0.8,  True, False)

executable_name = os.path.basename(executable)

answer1 = input("Do you want to launch the executable? [y/n]: ")
ask_execute = (answer1=="y")


if executable_name == "advection_ALL":
    answer2 = input(explanations + "Which tests do you want to display ? [int between 0 and 47, separated by \',\' or  \'-\']: ")
    list_string = answer2.split(",")

answer3 = input("Do you want to remove the folder containing all the data of each test case after running ? [y/n]: ")
ask_remove = (answer3=="y")

if ask_execute :
    # Execute the test file given as input in the command ----------
    out = execute(executable, yaml_parameters, verbose)


# Display the curves --------------------------------------------
if executable_name == "advection_ALL":
    if "-" not in answer2 :
        selected_test = [int(el) for el in list_string]
    else :
        selected_test = []
        for el in list_string:
            if "-" not in el:
                selected_test += [int(el)]
            else :
                [a, b] = el.split("-")
                selected_test += list(range(int(a), int(b) + 1))
    possible_output_folder_names = [d for d in os.listdir() if d.endswith('_output') and d.count('-')==2]
    output_folder_names = [possible_output_folder_names[s] for s in selected_test]
    keys = [folder.replace('_output','').split('-') for folder in output_folder_names]
    fct_names = [f'{simulation.capitalize()} with {method.replace("_"," ")} on {mapping.replace("_", " ")}'
                 for mapping, method, simulation in keys]
else :
    mapping, method, domain, simulation, name = get_simulation_config(executable)
    output_folder_names = [f'{mapping.replace(" ","_")}_{domain}-{method.replace(" ","_")}-{simulation}_output']
    fct_names = [name]
    selected_test = [-1]

rmin = yaml_parameters['SplineMesh']['r_min']
rmax = yaml_parameters['SplineMesh']['r_max']
Nr = yaml_parameters['SplineMesh']['r_ncells']
Nt = yaml_parameters['SplineMesh']['theta_ncells']
dt0 = yaml_parameters['Time']['time_step']
T = yaml_parameters['Time']['final_time']
details1 = f"_{Nr}x{Nt}_[{rmin}_{rmax}]"
for i, (name, folder) in enumerate(zip(fct_names, output_folder_names)):
    dt = dt0
    if 'euler' in folder:
        dt *= 0.1
    iter_nb = int(T/dt)
    details2 = f"\n $NrxNt$ = {Nr}x{Nt}; [$rmin$,$rmax$] = [{rmin},{rmax}]; dt = {dt}"
    animate(iter_nb, name.replace(" ", "_")+details1, name+details2, folder)

    if executable_name == "advection_ALL":
        print(f"The animation of the {selected_test} test case have been saved in ./animations/ folder.")
    else:
        print("The animation of the test case have been saved in ./animations/ folder.")




if ask_remove:
    delete_folders = possible_output_folder_names if executable_name == "advection_ALL" else output_folder_names
    for folder in delete_folders:
        shutil.rmtree(folder)
        print(f"The {folder} folder has been removed.")


