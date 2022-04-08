'''
file : POPEread.py
date : 2022/02/17

Reading of the POPE result files
'''

# pylint: disable=invalid-name

import numpy as np

import HDF5utils as H5ut

# pylint: disable=invalid-name

# ===============================================
#  Select the results required for POPE analysis
# -----------------------------------------------
def Read_POPE_results(H, ifirst_diag, iend_diag):
    '''
    Select the results required for POPE analysis
    '''

    POPE_read = H5ut.HDF5Group()

    # ---> Reading of the mesh
    POPE_read.xgrid = H.MeshX
    POPE_read.vgrid = H.MeshVx
    POPE_read.nx = H.Nx
    POPE_read.nv = H.Nvx

    # ---> Reading of the mass of ions and electrons
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
    POPE_read.sqrt_me_on_mi = np.sqrt(me/mi)

    # ---> Reading of the mask
    if hasattr(H, 'mask_rho'):
        POPE_read.mask_rho = H.mask_rho
    else:
        POPE_read.mask_rho = np.ones((1, POPE_read.nx))

    # ---> Reading of the time evolution of the potential and
    # --->  the distribution functions
    POPE_read.timegrid = H.time_saved[ifirst_diag:iend_diag]
    POPE_read.nt = len(POPE_read.timegrid)

    POPE_read.Phi = H.electrostatic_potential[ifirst_diag:iend_diag, :]
    POPE_read.felec = H.fdistribu[ifirst_diag:iend_diag, ielec, :, :]
    nb_kin_species = np.shape(H.fdistribu)[1]
    if nb_kin_species > 1:
        POPE_read.fion = H.fdistribu[ifirst_diag:iend_diag, ion, :, :]

    return POPE_read

# end def Read_POPE_results
