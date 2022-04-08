'''
POPE analysis for Vlasov-Poisson
   Author: V. Grandgirard
'''

# pylint: disable=invalid-name

import os

import numpy as np
import scipy.linalg as la

import HDF5utils as H5ut

# ------------------------------------------------------
#  Solve Least Square problem
#   y = Xb + e where
#    . X is a (N,k) matrix with N>>k
#    . e is the vector of residuals and
#    . b is the vector of estimates
#   Then the least square estimator is
#     b = (X'X)^{-1}X'y where X' is the transpose
#   matrix of X and providing that the inverse of
#   X'X exist (i.e X is of rank k).
#
#   i.e b = X^{+}y where X^{+} = (X'X)^{-1}X'
#   is the pseudo-inverse matrix of X
# ------------------------------------------------------
def solve_LeastSquare_(MLS, op_RHS_slice, i1, i2, *args):
    '''
    Solve Least Square problem
    y = Xb + e where
    . X is a (N,k) matrix with N>>k
    . e is the vector of residuals and
    . b is the vector of estimates
    Then the least square estimator is
    b = (X'X)^{-1}X'y where X' is the transpose
    matrix of X and providing that the inverse of
    X'X exist (i.e X is of rank k).

    i.e b = X^{+}y where X^{+} = (X'X)^{-1}X'
    is the pseudo-inverse matrix of X
    '''
    residual = 0.
    MLS_rank = np.linalg.matrix_rank(MLS)
    if MLS_rank == MLS.shape[1]:
        MLS_pseudoinv = np.linalg.pinv(MLS)
        xLS = np.matmul(MLS_pseudoinv, op_RHS_slice)
        for i, coeff_LHS in enumerate(args):
            coeff_LHS[i1, i2] = xLS[i]

        residual = np.mean(op_RHS_slice - np.matmul(MLS, xLS))

    return residual

# end def solve_LeastSquare_


# --------------------------------------------------
#  Computation of POPE coefficients for an equation
#   with 2 operators
# --------------------------------------------------
def compute_POPEcoeff_2operators(op_LHS_1, op_LHS_2, op_RHS, msg=''):
    '''
    Computation of POPE coefficients for an equation
    with 2 operators
    '''

    [nt, nx, nv] = np.shape(op_RHS)

    # --- Initialisation of the Projection matrix ---
    nb_op = 2
    M_pope = np.zeros((nb_op, nb_op))
    RHS_pope = np.zeros(nb_op)

    # --- Initialisation of the matrix for the projection determinant ---
    Det_pope = np.zeros([nt, nx, nv])

    # --- Initizalisation of the POPE result matrix  ---
    coeff_LHS1 = np.zeros([nt, nx, nv])
    coeff_LHS2 = np.zeros([nt, nx, nv])
    coeff_LHS1_LS_x = np.zeros([nt, nv])
    coeff_LHS2_LS_x = np.zeros([nt, nv])
    coeff_LHS1_LS_v = np.zeros([nt, nx])
    coeff_LHS2_LS_v = np.zeros([nt, nx])

    for it in range(nt):
        print(msg, " it = ", it, " / ", nt-1)
        for ix in range(nx-1):
            for iv in range(nv):
                M_pope = np.array([[op_LHS_1[it, ix, iv], op_LHS_2[it, ix, iv]],
                                   [op_LHS_1[it, ix+1, iv], op_LHS_2[it, ix+1, iv]]])
                Det_pope[it, ix, iv] = la.det(M_pope)
                if abs(Det_pope[it, ix, iv]) >= 1e-14:
                    RHS_pope = np.array([op_RHS[it, ix, iv], op_RHS[it, ix+1, iv]])
                    x_pope = la.solve(M_pope, RHS_pope)
                    coeff_LHS1[it, ix, iv] = x_pope[0]
                    coeff_LHS2[it, ix, iv] = x_pope[1]

        # Least square computation in x direction
        for iv in range(nv):
            MLS = np.array([op_LHS_1[it, :, iv], op_LHS_2[it, :, iv]])
            solve_LeastSquare_(np.transpose(MLS), op_RHS[it, :, iv], it, iv,
                               coeff_LHS1_LS_x, coeff_LHS2_LS_x)

        # Least square computation in v direction
        if nv > 1:
            for ix in range(nx):
                MLS = np.array([op_LHS_1[it, ix, :], op_LHS_2[it, ix, :]])
                solve_LeastSquare_(np.transpose(MLS), op_RHS[it, ix, :],
                                   it, ix, coeff_LHS1_LS_v, coeff_LHS2_LS_v)

    if nv > 1:
        return [coeff_LHS1, coeff_LHS2,
                coeff_LHS1_LS_x, coeff_LHS2_LS_x,
                coeff_LHS1_LS_v, coeff_LHS2_LS_v]
    else:
        return [coeff_LHS1[:, :, 0], coeff_LHS2[:, :, 0],
                coeff_LHS1_LS_x[:, 0], coeff_LHS2_LS_x[:, 0]]

# end def compute_POPEcoeff_2operators


# --------------------------------------------------
#  Computation of POPE coefficients for an equation
#   with 4 operators
# --------------------------------------------------
def compute_POPEcoeff_4operators(op_LHS_1, op_LHS_2,
                                 op_LHS_3, op_LHS_4, op_RHS, msg=''):
    '''
    Computation of POPE coefficients for an equation
    with 4 operators
    '''
    [nt, nx, nv] = np.shape(op_RHS)

    # --- Initialisation of the Projection matrix ---
    nb_op = 4
    M_pope = np.zeros((nb_op, nb_op))
    RHS_pope = np.zeros(nb_op)

    # --- Initialisation of the matrix for the projection determinant ---
    Det_pope = np.zeros([nt, nx, nv])

    # --- Initizalisation of the POPE result matrix  ---
    coeff_LHS1 = np.zeros([nt, nx, nv])
    coeff_LHS2 = np.zeros([nt, nx, nv])
    coeff_LHS3 = np.zeros([nt, nx, nv])
    coeff_LHS4 = np.zeros([nt, nx, nv])
    coeff_LHS1_LS_x = np.zeros([nt, nv])
    coeff_LHS2_LS_x = np.zeros([nt, nv])
    coeff_LHS3_LS_x = np.zeros([nt, nv])
    coeff_LHS4_LS_x = np.zeros([nt, nv])
    coeff_LHS1_LS_v = np.zeros([nt, nx])
    coeff_LHS2_LS_v = np.zeros([nt, nx])
    coeff_LHS3_LS_v = np.zeros([nt, nx])
    coeff_LHS4_LS_v = np.zeros([nt, nx])

    for it in range(nt):
        print(msg, " op it = ", it, " / ", nt-1)
        for iv in range(nv-1):
            for ix in range(nx-1):
                M_pope = np.array(
                    [[op_LHS_1[it, ix, iv], op_LHS_2[it, ix, iv],
                      op_LHS_3[it, ix, iv], op_LHS_4[it, ix, iv]],
                     [op_LHS_1[it, ix+1, iv], op_LHS_2[it, ix+1, iv],
                      op_LHS_3[it, ix+1, iv], op_LHS_4[it, ix+1, iv]],
                     [op_LHS_1[it, ix, iv+1], op_LHS_2[it, ix, iv+1],
                      op_LHS_3[it, ix, iv+1], op_LHS_4[it, ix, iv+1]],
                     [op_LHS_1[it, ix+1, iv+1], op_LHS_2[it, ix+1, iv+1],
                      op_LHS_3[it, ix+1, iv+1], op_LHS_4[it, ix+1, iv+1]]])
                Det_pope[it, ix, iv] = la.det(M_pope)
                if abs(Det_pope[it, ix, iv]) >= 1e-14:
                    RHS_pope = np.array(
                        [op_RHS[it, ix, iv],
                         op_RHS[it, ix+1, iv],
                         op_RHS[it, ix, iv+1],
                         op_RHS[it, ix+1, iv+1]])
                    x_pope = la.solve(M_pope, RHS_pope)
                    coeff_LHS1[it, ix, iv] = x_pope[0]
                    coeff_LHS2[it, ix, iv] = x_pope[1]
                    coeff_LHS3[it, ix, iv] = x_pope[2]
                    coeff_LHS4[it, ix, iv] = x_pope[3]

        # Least square computation in x direction
        for iv in range(nv):
            MLS = np.array(
                [op_LHS_1[it, :, iv],
                 op_LHS_2[it, :, iv],
                 op_LHS_3[it, :, iv],
                 op_LHS_4[it, :, iv]])
            solve_LeastSquare_(np.transpose(MLS), op_RHS[it, :, iv], it, iv,
                               coeff_LHS1_LS_x, coeff_LHS2_LS_x,
                               coeff_LHS3_LS_x, coeff_LHS4_LS_x)

        # Least square computation in v direction
        for ix in range(nx):
            MLS = np.array(
                [op_LHS_1[it, ix, :],
                 op_LHS_2[it, ix, :],
                 op_LHS_3[it, ix, :],
                 op_LHS_4[it, ix, :]])
            solve_LeastSquare_(np.transpose(MLS), op_RHS[it, ix, :], it, ix,
                               coeff_LHS1_LS_v, coeff_LHS2_LS_v,
                               coeff_LHS3_LS_v, coeff_LHS4_LS_v)

    return [coeff_LHS1, coeff_LHS2, coeff_LHS3, coeff_LHS4,
            coeff_LHS1_LS_x, coeff_LHS2_LS_x, coeff_LHS3_LS_x, coeff_LHS4_LS_x,
            coeff_LHS1_LS_v, coeff_LHS2_LS_v, coeff_LHS3_LS_v, coeff_LHS4_LS_v]

# end def compute_POPEcoeff_4operators


# --------------------------------------------------
#  Computation of the POPE coefficients
#   for Vlasov-Poisson problem without RHS
# --------------------------------------------------
def Compute_POPEcoeff_forPoisson(
        POPE_read,
        dic_POPEop_Poisson):
    '''
    Computation of the POPE coefficients for Poisson equation
      -d2Phi/dx2 = rho
        where rho(x) = Mp(x)*(qi*intdv_fion - qe*intdv_felec)
        with Mp(x) a mask to delimit plasma region
    '''

    mask_rho = POPE_read.mask_rho

    # --- Computation of qe and qi ---
    d2xPhi = dic_POPEop_Poisson['d2xPhi']
    ni = dic_POPEop_Poisson['ni']
    ne = dic_POPEop_Poisson['ne']

    Mp = mask_rho
    # --> Poisson : op_LHS_1 + op_LHS_2 = op_RHS
    op_LHS_1 = Mp*ni
    op_LHS_2 = -Mp*ne
    op_RHS = -d2xPhi
    [qi, qe, qi_LS, qe_LS] = compute_POPEcoeff_2operators(
        op_LHS_1[:, :, None], op_LHS_2[:, :, None], op_RHS[:, :, None], "Poisson")

    dic_POPEcoeff_Poisson = {
        'xgrid': dic_POPEop_Poisson['xgrid'],
        'timegrid': dic_POPEop_Poisson['timegrid'],
        'qe': qe,
        'qi': qi,
        'qe_LS': qe_LS,
        'qi_LS': qi_LS}

    return dic_POPEcoeff_Poisson

# end def Compute_POPEcoeff_forPoisson


# --------------------------------------------------
#  Get (load if already exist or compute)
#   the POPE coefficients for the Poisson equation
# --------------------------------------------------
def Get_POPEcoeff_forPoisson(
        POPE_read,
        dic_POPEop_Poisson,
        POPE_h5filename):
    '''
    Get (load if already exist or compute)
     the POPE coefficients for the Poisson equation
    '''

    # ---> Check if the POPE coefficients for Poisson are already saved
    h5group_POPEcoeff_Poisson = 'POPEcoeff_Poisson/'
    if os.path.exists(POPE_h5filename):
        dic_POPEcoeff_Poiss_exist = H5ut.group_exist_in_hdf5(
            POPE_h5filename, h5group_POPEcoeff_Poisson)
        if dic_POPEcoeff_Poiss_exist:
            print('==> Loading of POPE coefficients for Poisson in '+POPE_h5filename)
            dic_POPEcoeff_Poisson = H5ut.load_dict_from_hdf5(
                POPE_h5filename, h5group_POPEcoeff_Poisson)
            return dic_POPEcoeff_Poisson

    print('==> Computation of POPE coefficients for Poisson')
    # ---> Computation of PoPe coefficients for Poisson
    dic_POPEcoeff_Poisson = Compute_POPEcoeff_forPoisson(
        POPE_read, dic_POPEop_Poisson)

    print('==> Saving of POPE coefficients for Poisson in '+POPE_h5filename)
    H5ut.save_dict_to_hdf5(
        dic_POPEcoeff_Poisson, h5group_POPEcoeff_Poisson, POPE_h5filename)

    return dic_POPEcoeff_Poisson

# end def Get_POPEcoeff_forPoisson


# --------------------------------------------------
#  Computation of the POPE coefficients
#   for Vlasov equations without RHS
# --------------------------------------------------
def Compute_POPEcoeff_forVlasov_noRHS(
        POPE_read, dic_POPEop_Vlasov):
    '''
     Computation of the POPE coefficients for the Vlasov equation
      dfelec/dt = - ae*v*dfelec/dx - be*dPhi/dx*dfelec/dv
      dfion/dt  = - ai*sqrt_me_on_mi*v*dfion/dx
                  + bi*sqrt_me_on_mi*dPhi/dx*dfion/dv
    '''

    sqrt_me_on_mi = POPE_read.sqrt_me_on_mi

    # --- Computation of ae,be, ai, bi ---
    # ------- Vlasov electrons ------
    dxPhi = dic_POPEop_Vlasov['dxPhi']
    v_dxfelec = dic_POPEop_Vlasov['v_dxfelec']
    dvfelec = dic_POPEop_Vlasov['dvfelec']
    dtfelec = dic_POPEop_Vlasov['dtfelec']

    op_LHS_1 = -v_dxfelec
    op_LHS_2 = -dxPhi[:, :, None]*dvfelec
    op_RHS = dtfelec

    [ae, be, ae_LS_x, be_LS_x, ae_LS_v, be_LS_v] = compute_POPEcoeff_2operators(
        op_LHS_1, op_LHS_2, op_RHS, "Vlasov elec")

    # --------- Vlasov ions ---------
    if hasattr(POPE_read, 'fion'):
        dxPhi = dic_POPEop_Vlasov['dxPhi']
        v_dxfion = dic_POPEop_Vlasov['v_dxfion']
        dvfion = dic_POPEop_Vlasov['dvfion']
        dtfion = dic_POPEop_Vlasov['dtfion']

        op_LHS_1 = -sqrt_me_on_mi*v_dxfion
        op_LHS_2 = sqrt_me_on_mi*dxPhi[:, :, None]*dvfion
        op_RHS = dtfion

        [ai, bi, ai_LS_x, bi_LS_x, ai_LS_v, bi_LS_v] = compute_POPEcoeff_2operators(
            op_LHS_1, op_LHS_2, op_RHS, "Vlasov ion")

    dic_POPEcoeff_Vlasov = {
        'timegrid': dic_POPEop_Vlasov['timegrid'],
        'xgrid': dic_POPEop_Vlasov['xgrid'],
        'vgrid': dic_POPEop_Vlasov['vgrid'],
        'ae': ae,
        'be': be,
        'ae_LS_x': ae_LS_x,
        'be_LS_x': be_LS_x,
        'ae_LS_v': ae_LS_v,
        'be_LS_v': be_LS_v}
    if hasattr(POPE_read, 'fion'):
        dic_POPEcoeff_Vlasov['ai'] = ai
        dic_POPEcoeff_Vlasov['bi'] = bi
        dic_POPEcoeff_Vlasov['ai_LS_x'] = ai_LS_x
        dic_POPEcoeff_Vlasov['bi_LS_x'] = bi_LS_x
        dic_POPEcoeff_Vlasov['ai_LS_v'] = ai_LS_v
        dic_POPEcoeff_Vlasov['bi_LS_v'] = bi_LS_v

    return dic_POPEcoeff_Vlasov

# end def Compute_POPEcoeff_forVlasov_noRHS


# --------------------------------------------------
#  Get (load if already exist or compute)
#   the POPE coefficients for the Vlasov equation
#   including RHS if existing
# --------------------------------------------------
def Get_POPEcoeff_forVlasov(
        POPE_read,
        dic_POPEop_Vlasov,
        POPE_h5filename):
    '''
    Get (load if already exist or compute)
     the POPE coefficients for the Vlasov equation
     including RHS if existing
    '''

    # ---> Check if the POPE coefficients are already saved
    h5group_POPEcoeff_Vlasov = 'POPEcoeff_Vlasov/'
    if os.path.exists(POPE_h5filename):
        dic_POPEcoeff_Vlasov_exist = H5ut.group_exist_in_hdf5(
            POPE_h5filename, h5group_POPEcoeff_Vlasov)
        if dic_POPEcoeff_Vlasov_exist:
            print('==> Loading of POPE coefficients for Vlasov in '+POPE_h5filename)
            dic_POPEcoeff_Vlasov = H5ut.load_dict_from_hdf5(
                POPE_h5filename, h5group_POPEcoeff_Vlasov)
            return dic_POPEcoeff_Vlasov

    print('==> Computation of POPE coefficients for Vlasov')
    dic_POPEcoeff_Vlasov = Compute_POPEcoeff_forVlasov_noRHS(
        POPE_read, dic_POPEop_Vlasov)

    print('==> Saving of POPE coefficients for Vlasov in '+POPE_h5filename)
    H5ut.save_dict_to_hdf5(
        dic_POPEcoeff_Vlasov, h5group_POPEcoeff_Vlasov, POPE_h5filename)

    return dic_POPEcoeff_Vlasov

# end def Get_POPEcoeff_forVlasov


# --------------------------------------------------
# --> Statistic analysis on POPE coefficients
#      Build the Gaussian
# --------------------------------------------------
def ComputeGaussian_POPE_coeff(dic_POPE_coeff):
    ''''
    POPE error statistic analysis
    '''

    for key_name in dic_POPE_coeff:
        if key_name not in ['xgrid', 'vgrid', 'timegrid']:
            key_val = dic_POPE_coeff[key_name]
            key_val_nonzero = key_val[key_val != 0.]
            if len(key_val_nonzero) != 0:
                key_mean = np.mean(key_val_nonzero)
                key_RMS = np.sqrt(np.mean((key_val_nonzero-key_mean) *
                                          (key_val_nonzero-key_mean)))
                print('{}: Mean = {}, RMS = {}'.format(key_name, key_mean, key_RMS))
            else:
                print('{}: Mean = NAN, RMS = NAN'.format(key_name))

# end def ComputeGaussian_POPE_coeff
