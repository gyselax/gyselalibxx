# --------------------------------------------------
#  iPOPE analysis for Vlasov-Poisson
#    Author: V. Grandgirard
# --------------------------------------------------
'''
iPOPE analysis for Vlasov-Poisson
'''

# pylint: disable=invalid-name

import os

import numpy as np

import HDF5utils as H5ut
import MATHutils as MATHut

# --------------------------------------------------
#  Computation of the POPE error and residu
#   for Poisson equation
# --------------------------------------------------
def Compute_iPOPEcoeff_forPoisson(
        POPE_read,
        dic_POPEop_Poisson):
    '''
    Computation of the IPOPE error and residu
     -d2Phi/dx2 = rho
        where rho(x) = Mp(x)*(qi*intdv_fion - qe*intdv_felec)
        with Mp(x) a mask to delimit plasma region
    '''
    mask_rho = POPE_read.mask_rho

    timegrid = dic_POPEop_Poisson['timegrid']
    xgrid = dic_POPEop_Poisson['xgrid']
    vgrid = dic_POPEop_Poisson['vgrid']
    nt = len(timegrid)
    nx = len(xgrid)

    weight_q = np.ones([nx])
    rhs_q = np.zeros([nx])
    iPOPE_error_q = np.zeros([nx])
    Coef_q_num = np.zeros([nx])
    Coef_q_denom = np.full([nx], 1.e-16)
    iPOPE_index_q = np.zeros([nx])

    d2xPhi = dic_POPEop_Poisson['d2xPhi']
    ni = dic_POPEop_Poisson['ni']
    ne = dic_POPEop_Poisson['ne']

    Mp = mask_rho
    for it in range(nt):
        print("iPOPE for Poisson at time " + str(it) + "--> " + str(timegrid[it]))
        rhs_q_i = Mp[:]*(ni[it, :] - ne[it, :])
        rhs_q[:] = rhs_q[:] + rhs_q_i
        error_q_i = -d2xPhi[it, :] - rhs_q_i

        # --> 0: weight=1. (no more used 1:weight=|rhs|, 2:weight=1/|rhs|)
        weight_q_i = 1.
        weight_q[:] = weight_q[:] + weight_q_i
        iPOPE_error_q[:] = iPOPE_error_q[:] + error_q_i*weight_q_i

        Coef_q_num[:] = Coef_q_num[:] + error_q_i*rhs_q_i*weight_q_i
        Coef_q_denom[:] = Coef_q_denom[:] + rhs_q_i*rhs_q_i*weight_q_i

    iPOPE_Coef_q = Coef_q_num / Coef_q_denom
    iPOPE_index_q = -np.log10(MATHut.replaceZeros(np.abs(iPOPE_Coef_q)))
    iPOPE_error_q = iPOPE_error_q / weight_q
    iPOPE_residu_q = (iPOPE_error_q - iPOPE_Coef_q*rhs_q/weight_q)

    dic_iPOPEcoeff_Poisson = {
        'xgrid': xgrid,
        'vgrid': vgrid,
        'timegrid': timegrid,
        'iPOPE_index_q': iPOPE_index_q,
        'iPOPE_error_q': iPOPE_error_q,
        'iPOPE_residu_q': iPOPE_residu_q}

    return dic_iPOPEcoeff_Poisson

# end def Compute_iPOPEcoeff_forPoisson


# ---------------------------------------------------------
#  Get (load if already exist or compute)
#   the iPOPE error and residu for the Poisson equation
# ---------------------------------------------------------
def Get_iPOPEcoeff_forPoisson(
        POPE_read,
        dic_POPEop_Poisson,
        POPE_h5filename):
    '''
    Get (load if already exist or compute)
     the iPOPE error and residu for the Poisson equation
    '''

    # ---> Check if the iPOPE coefficients are already saved
    h5group_iPOPEcoeff_Poisson = 'iPOPEcoeff_Poisson/'
    if os.path.exists(POPE_h5filename):
        dic_iPOPEcoeff_Poisson_exist = H5ut.group_exist_in_hdf5(
            POPE_h5filename, h5group_iPOPEcoeff_Poisson)
        if dic_iPOPEcoeff_Poisson_exist:
            print('==> Loading of iPOPE coefficients for Poisson in ' +
                  POPE_h5filename)
            dic_iPOPEcoeff_Poisson = H5ut.load_dict_from_hdf5(
                POPE_h5filename, h5group_iPOPEcoeff_Poisson)
            return dic_iPOPEcoeff_Poisson

    print('==> Saving of iPOPE coefficients for Poisson in ' + POPE_h5filename)
    dic_iPOPEcoeff_Poisson = Compute_iPOPEcoeff_forPoisson(
        POPE_read, dic_POPEop_Poisson)
    H5ut.save_dict_to_hdf5(dic_iPOPEcoeff_Poisson,
                           h5group_iPOPEcoeff_Poisson, POPE_h5filename)

    return dic_iPOPEcoeff_Poisson

# end Get_iPOPEcoeff_forPoisson


# --------------------------------------------------
#  Computation of the iPOPE error and residu
#   for the Vlasov equation
# --------------------------------------------------
def Compute_iPOPEcoeff_forVlasov(
        POPE_read, dic_POPEop_Vlasov):
    '''
    Computation of the IPOPE error and residu
      C = sum_j [df/dt-V_j
      E = -dPhi/dx
      dfelec/dt = -v*dfelec/dx + E*dfelec/dv
      dfion/dt  = - sqrt(me/mi)*(v*dfion/dx+E*dfion/dv)
    '''
    sqrt_me_on_mi = POPE_read.sqrt_me_on_mi

    timegrid = dic_POPEop_Vlasov['timegrid']
    xgrid = dic_POPEop_Vlasov['xgrid']
    vgrid = dic_POPEop_Vlasov['vgrid']
    nt = len(timegrid)
    nx = len(xgrid)
    nv = len(vgrid)

    dxPhi = dic_POPEop_Vlasov['dxPhi']

    weight_e = np.ones([nx, nv])
    rhs_e = np.zeros([nx, nv])
    iPOPE_error_e = np.zeros([nx, nv])
    Coef_e_num = np.zeros([nx, nv])
    Coef_e_denom = np.full([nx, nv], 1.e-16)
    iPOPE_index_e = np.zeros([nx, nv])

    # seuil_min = np.full([nx, nv], 1.e-16)

    # ---> for elec
    v_dxfelec = dic_POPEop_Vlasov['v_dxfelec']
    dvfelec = dic_POPEop_Vlasov['dvfelec']
    dtfelec = dic_POPEop_Vlasov['dtfelec']

    for it in range(nt):
        print("iPOPE for Vlasov for electrons at time " +
              str(it) + "--> " + str(timegrid[it]))
        for ix in range(1, nx):
            rhs_e_ij = - v_dxfelec[it, ix, :] - \
                dxPhi[it, ix, None]*dvfelec[it, ix, :]
            rhs_e[ix, :] = rhs_e[ix, :] + rhs_e_ij
            error_e_ij = dtfelec[it, ix, :] - rhs_e_ij
            # --> 0: weight=1. (no more used 1:weight=|rhs|, 2:weight=1/|rhs|)
            weight_e_ij = 1.
            weight_e[ix, :] = weight_e[ix, :] + weight_e_ij
            iPOPE_error_e[ix, :] = iPOPE_error_e[ix, :] + error_e_ij*weight_e_ij
            Coef_e_num[ix, :] = Coef_e_num[ix, :] + \
                error_e_ij*rhs_e_ij*weight_e_ij
            Coef_e_denom[ix, :] = Coef_e_denom[ix, :] + \
                rhs_e_ij*rhs_e_ij*weight_e_ij
        # end for ix
    # end for it
    iPOPE_Coef_e = Coef_e_num / Coef_e_denom
    iPOPE_index_e = -np.log10(MATHut.replaceZeros(np.abs(iPOPE_Coef_e)))
    iPOPE_error_e = iPOPE_error_e / weight_e
    iPOPE_residu_e = (iPOPE_error_e - iPOPE_Coef_e*rhs_e/weight_e)

    # ---> for ion
    if hasattr(POPE_read, 'fion'):
        weight_i = np.ones([nx, nv])
        rhs_i = np.zeros([nx, nv])
        iPOPE_error_i = np.zeros([nx, nv])
        Coef_i_num = np.zeros([nx, nv])
        Coef_i_denom = np.full([nx, nv], 1.e-16)
        iPOPE_index_i = np.zeros([nx, nv])

        v_dxfion = dic_POPEop_Vlasov['v_dxfion']
        dvfion = dic_POPEop_Vlasov['dvfion']
        dtfion = dic_POPEop_Vlasov['dtfion']

        for it in range(nt):
            print("iPOPE for Vlasov at time " + str(it) + "--> " + str(timegrid[it]))
            for ix in range(1, nx):
                rhs_i_ij = - sqrt_me_on_mi * (
                    v_dxfion[it, ix, :] - dxPhi[it, ix, None]*dvfion[it, ix, :])
                rhs_i[ix, :] = rhs_i[ix, :] + rhs_i_ij
                error_i_ij = dtfion[it, ix, :] - rhs_i_ij
                # --> 0: weight=1. (no more used 1:weight=|rhs|, 2:weight=1/|rhs|)
                weight_i_ij = 1.
                weight_i[ix, :] = weight_i[ix, :] + weight_i_ij
                iPOPE_error_i[ix, :] = iPOPE_error_i[ix, :] + \
                    error_i_ij*weight_i_ij
                Coef_i_num[ix, :] = Coef_i_num[ix, :] + \
                    error_i_ij*rhs_i_ij*weight_i_ij
                Coef_i_denom[ix, :] = Coef_i_denom[ix, :] + \
                    rhs_i_ij*rhs_i_ij*weight_i_ij
            # end for ix
        # end for it
        iPOPE_Coef_i = Coef_i_num / Coef_i_denom
        iPOPE_index_i = -np.log10(MATHut.replaceZeros(np.abs(iPOPE_Coef_i)))
        iPOPE_error_i = iPOPE_error_i / weight_i
        iPOPE_residu_i = (iPOPE_error_i - iPOPE_Coef_i*rhs_i/weight_i)

    dic_iPOPEcoeff_Vlasov = {
        'xgrid': xgrid,
        'vgrid': vgrid,
        'timegrid': timegrid,
        'iPOPE_index_e': iPOPE_index_e,
        'iPOPE_error_e': iPOPE_error_e,
        'iPOPE_residu_e': iPOPE_residu_e}
    if hasattr(POPE_read, 'fion'):
        dic_iPOPEcoeff_Vlasov['iPOPE_index_i'] = iPOPE_index_i
        dic_iPOPEcoeff_Vlasov['iPOPE_error_i'] = iPOPE_error_i
        dic_iPOPEcoeff_Vlasov['iPOPE_residu_i'] = iPOPE_residu_i

    return dic_iPOPEcoeff_Vlasov

# end Compute_iPOPEcoeff_forVlasov


# ----------------------------------------------------
#  Get (load if already exist or compute)
#   the iPOPE error and residu for the Vlasov equation
# ----------------------------------------------------
def Get_iPOPEcoeff_forVlasov(
        POPE_read,
        dic_POPEop_Vlasov,
        POPE_h5filename):
    '''
    Get (load if already exist or compute)
     the iPOPE error and residu for the Vlasov equation
    '''

    # ---> Check if the iPOPE coefficients are already saved
    h5group_iPOPEcoeff_Vlasov = 'iPOPEcoeff_Vlasov/'
    if os.path.exists(POPE_h5filename):
        dic_iPOPEcoeff_Vlasov_exist = H5ut.group_exist_in_hdf5(
            POPE_h5filename, h5group_iPOPEcoeff_Vlasov)
        if dic_iPOPEcoeff_Vlasov_exist:
            print('==> Loading of iPOPE coefficients in ' + POPE_h5filename)
            dic_iPOPEcoeff_Vlasov = H5ut.load_dict_from_hdf5(
                POPE_h5filename, h5group_iPOPEcoeff_Vlasov)
            return dic_iPOPEcoeff_Vlasov

    # ---> Computation of iPoPe coefficients for Vlasov
    print('==> Saving of iPOPE coefficients for Vlasov in ' + POPE_h5filename)
    dic_iPOPEcoeff_Vlasov = Compute_iPOPEcoeff_forVlasov(
        POPE_read, dic_POPEop_Vlasov)
    H5ut.save_dict_to_hdf5(dic_iPOPEcoeff_Vlasov,
                           h5group_iPOPEcoeff_Vlasov, POPE_h5filename)

    return dic_iPOPEcoeff_Vlasov

# end Get_iPOPEcoeff_forVlasov
