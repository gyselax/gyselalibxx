"""
Functions used to read the PINN parameters
and the results coming from VOICEXX
"""

# pylint: disable=E1101
# pylint: disable=invalid-name
# pylint: disable=import-error

import numpy as np
import yaml

import HDF5utils as H5ut
import MATHutils as MATHut
import READutils as READut

#-----------------------------------------------
# Read the PINN parameter from the YAML input file
#-----------------------------------------------
def Read_PINN_parameters(YAML_inputFile):
    '''
    Read the PINN parameter from the YAML input file
    '''

    with open(YAML_inputFile, "r") as stream:
        input_data = yaml.safe_load(stream)

    for key, data in input_data.items():
        if isinstance(data, str):
            if data.lower() == 'none':
                input_data[key] = None

    return input_data

#end def Read_PINN_parameters


#-----------------------------------------------
# Select the results from VOICEXX required for
# PINN analysis
#-----------------------------------------------
def Read_PINN_fromVOICEXX(VOICEXX_dir):
    '''
    Select the results from VOICEXX required for PINN analysis
    '''
    H = READut.Read_VOICEXX_results(VOICEXX_dir)

    PINN_read = H5ut.HDF5Group()

    #---> Reading of the mesh
    PINN_read.xgrid = H.MeshX
    PINN_read.vgrid = H.MeshVx
    PINN_read.nx = H.Nx
    PINN_read.nv = H.Nvx

    #---> Reading of the mass of ions and electrons
    fdistribu_charges = H.fdistribu_charges
    fdistribu_masses = H.fdistribu_masses
    [ielec_find] = np.where(fdistribu_charges == -1)
    ielec = ielec_find[0]
    [ion_find] = np.where(fdistribu_charges > 0)
    if np.size(ion_find) == 0:
        me = fdistribu_masses
        mi = 1.
    else:
        ion = ion_find[0]
        me = fdistribu_masses[ielec]
        mi = fdistribu_masses[ion]
    PINN_read.sqrt_me_on_mi = np.sqrt(me/mi)

    #---> Reading of the time evolution of the distribution functions
    PINN_read.timegrid = H.time_saved
    PINN_read.nt = len(PINN_read.timegrid)

    PINN_read.felec = H.fdistribu[:, ielec, :, :]
    nb_kin_species = np.shape(H.fdistribu)[1]
    if nb_kin_species > 1:
        PINN_read.fion = H.fdistribu[:, ion, :, :]

    #---> Reading of the time evolution of the electrostatic potential
    PINN_read.Phi = H.electrostatic_potential

    #---> Computation of the electric field
    #VG# [TODO] ATTENTION only valid for equidistant mesh
    dx = PINN_read.xgrid[1]-PINN_read.xgrid[0]
    PINN_read.Efield = np.zeros([PINN_read.nt, PINN_read.nx], dtype=float)
    for it in range(len(PINN_read.timegrid)):
        PINN_read.Efield[it, :] = -MATHut.Deriv_nthorder(PINN_read.Phi[it, :], dx, periodic=1)

    return PINN_read

#end def Read_PINN_fromVOICEXX
