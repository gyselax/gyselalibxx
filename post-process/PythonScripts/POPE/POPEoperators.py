'''
Computation of the POPE operators
'''

# pylint: disable=invalid-name

import os

import numpy as np

import HDF5utils as H5ut
import MATHutils as MATHut

# --------------------------------------------------
#  Computation of the operators appearing in
#      Poisson equation
#     Author: V. Grandgirard
# --------------------------------------------------
def Compute_POPEop_forPoisson(
        FD_order, # Finite Difference order
        POPE_read):
    '''
    Computation of Poisson operators
      (intdv_fion), (intdv_felec), (d2Phi/dx2)
      the x derivative are performed in the fourier space
    '''

    timegrid = POPE_read.timegrid
    xgrid = POPE_read.xgrid
    vgrid = POPE_read.vgrid
    nt = len(timegrid)
    nx = len(POPE_read.xgrid)
    Phi = POPE_read.Phi
    felec = POPE_read.felec
    if hasattr(POPE_read, 'fion'):
        fion = POPE_read.fion

    nb_ghost = int(FD_order/2)+1
    nt_op = nt-2*nb_ghost

    timegrid_P = np.zeros([nt_op])
    d2xPhi = np.zeros([nt_op, nx])

    for it_op in range(nt_op):
        it = it_op + nb_ghost
        print("Compute Poisson operators at it " + str(it_op))
        timegrid_P[it_op] = timegrid[it]
        [FFTPhi, k] = MATHut.Fourier1D(Phi[it, :], xgrid)
        FFT2Phi = -k*k*FFTPhi
        [d2xPhi[it_op, :], k] = MATHut.iFourier1D(FFT2Phi, k)

    # Here the potential Phi is recomputed as the result of Poisson
    #  equation. Only felec and fion are used as output of VOICE
    Phi_POPE = np.zeros([nt_op, nx])
    d2xPhi_POPE = np.zeros([nt_op, nx])
    ne = MATHut.compute_intFdv(felec[nb_ghost:nt-nb_ghost, :, :], vgrid)

    if hasattr(POPE_read, 'fion'):
        ni = MATHut.compute_intFdv(fion[nb_ghost:nt-nb_ghost, :, :], vgrid)
    else:
        ni = np.ones(np.shape(ne))
    rho = (ni - ne)

    for it_op in range(nt_op):
        [FFTrho, k] = MATHut.Fourier1D(rho[it_op, :], xgrid)
        kmiddle = int(len(k)/2)
        k2 = k*k
        k2[kmiddle] = 1.   # Fourier mode = 0
        FFTrho = FFTrho/k2
        FFTrho[kmiddle] = 0.
        [Phi_POPE[it_op, :], k] = MATHut.iFourier1D(FFTrho, k)
        [FFTPhi_POPE, k] = MATHut.Fourier1D(Phi_POPE[it_op, :], xgrid)
        FFT2Phi_POPE = -k*k*FFTPhi_POPE
        [d2xPhi_POPE[it_op, :], _] = MATHut.iFourier1D(FFT2Phi_POPE, k)

    dic_POPEop_Poisson = {}
    dic_POPEop_Poisson['timegrid'] = timegrid_P
    dic_POPEop_Poisson['xgrid'] = xgrid
    dic_POPEop_Poisson['vgrid'] = vgrid
    dic_POPEop_Poisson['d2xPhi'] = d2xPhi
    dic_POPEop_Poisson['ni'] = ni
    dic_POPEop_Poisson['ne'] = ne
    dic_POPEop_Poisson['rho'] = rho
    dic_POPEop_Poisson['Phi_POPE'] = Phi_POPE
    dic_POPEop_Poisson['d2xPhi_POPE'] = d2xPhi_POPE

    return dic_POPEop_Poisson

# end def Compute_POPEop_forPoisson


# --------------------------------------------------
#  Get (load if already exist or compute)
#   the operators appearing in Poisson equation
# --------------------------------------------------
def Get_POPEop_forPoisson(
        FD_order, # Finite Difference order
        POPE_read,
        POPE_h5filename):
    '''
    Load the Poisson operators saved in POPE_h5filename
    '''

    h5group_Poisson = 'Poisson/'
    if os.path.exists(POPE_h5filename):
        dic_POPEop_Poisson_exist = H5ut.group_exist_in_hdf5(POPE_h5filename,
                                                            h5group_Poisson)
        if dic_POPEop_Poisson_exist:
            print('==> Loading of POPE operators for Poisson in '+POPE_h5filename)
            dic_POPEop_Poisson = H5ut.load_dict_from_hdf5(POPE_h5filename,
                                                          h5group_Poisson)
            return dic_POPEop_Poisson

    dic_POPEop_Poisson = Compute_POPEop_forPoisson(
        FD_order, # Finite Difference order
        POPE_read)

    print('==> Saving of POPE operators for Poisson in '+POPE_h5filename)
    H5ut.save_dict_to_hdf5(dic_POPEop_Poisson, h5group_Poisson, POPE_h5filename)

    return dic_POPEop_Poisson

# end def Get_POPEop_forPoisson


# --------------------------------------------------
#  Computation of the operators appearing in
#      Vlasov equations for Landau case
#     Author: V. Grandgirard
# --------------------------------------------------
def Compute_POPEop_forVlasovLHS(
        FD_order, # Finite Difference order
        POPE_read):
    '''
    Computation of Vlasov operators
     (dfelec/dt) , (v*dfelec/dx) , (dfelec/dv) ,
     (dfion/dt) , (v*dfion/dx) , (dfion/dv) ,
     (dPhi/dx)
     the x derivative are performed in the Fourier space
    '''
    timegrid = POPE_read.timegrid
    xgrid = POPE_read.xgrid
    vgrid = POPE_read.vgrid
    nt = len(timegrid)
    nx = len(xgrid)
    nv = len(vgrid)
    Phi = POPE_read.Phi
    felec = POPE_read.felec
    if hasattr(POPE_read, 'fion'):
        fion = POPE_read.fion
        fion_exist = True
    else:
        fion = np.zeros(np.shape(felec))
        fion_exist = False

    dv = abs(vgrid[1] - vgrid[0])    # ATTENTION only valid for uniform mesh
    dt = timegrid[1] - timegrid[0]

    nb_ghost = int(FD_order/2)+1
    nt_op = nt-2*nb_ghost
    timegrid_V = np.zeros([nt_op])

    # --------- electric potential ---------
    dxPhi = np.zeros([nt_op, nx])
    for it_op in range(nt_op):
        it = it_op + nb_ghost
        [FFTPhi, k] = MATHut.Fourier1D(Phi[it, :], xgrid)
        FFTPhi = 1.j*k*FFTPhi
        [dxPhi[it_op, :], _] = MATHut.iFourier1D(FFTPhi, k)

    # --------- Vlasov electrons ---------
    v_dxfelec = np.zeros([nt_op, nx, nv])
    dvfelec = np.zeros([nt_op, nx, nv])
    dtfelec = np.zeros([nt_op, nx, nv])

    for it_op in range(nt_op):
        it = it_op + nb_ghost
        print("Compute Vlasov operators for electrons at it " + str(it_op))
        timegrid_V[it_op] = timegrid[it]
        for iv in range(nv):
            [FFTxfelec, kfelec] = MATHut.Fourier1D(felec[it, :, iv], xgrid)
            FFTxfelec = 1.j*kfelec*FFTxfelec
            [v_dxfelec[it_op, :, iv], _] = MATHut.iFourier1D(FFTxfelec, kfelec)
            v_dxfelec[it_op, :, iv] = vgrid[iv] * v_dxfelec[it_op, :, iv]
        for ix in range(nx):
            dvfelec[it_op, ix, :] = MATHut.Deriv1(felec[it, ix, :], dv,
                                                  periodic=0, order=FD_order)
    for iv in range(nv):
        for ix in range(nx):
            dtfelec_tmp = MATHut.Deriv1(felec[:, ix, iv], dt,
                                        periodic=0, order=FD_order)
            dtfelec[:, ix, iv] = dtfelec_tmp[nb_ghost:-nb_ghost]

    # --------- Vlasov ions ---------
    if fion_exist:
        v_dxfion = np.zeros([nt_op, nx, nv])
        dvfion = np.zeros([nt_op, nx, nv])
        dtfion = np.zeros([nt_op, nx, nv])

        for it_op in range(nt_op):
            it = it_op + nb_ghost
            print("Compute Vlasov operators for ions at it " + str(it_op))
            timegrid_V[it_op] = timegrid[it]
            for iv in range(nv):
                [FFTxfion, kfion] = MATHut.Fourier1D(fion[it, :, iv], xgrid)
                FFTxfion = 1.j*kfion*FFTxfion
                [v_dxfion[it_op, :, iv], _] = MATHut.iFourier1D(FFTxfion, kfion)
                v_dxfion[it_op, :, iv] = vgrid[iv] * v_dxfion[it_op, :, iv]
            for ix in range(nx):
                dvfion[it_op, ix, :] = MATHut.Deriv1(fion[it, ix, :], dv,
                                                     periodic=0, order=FD_order)
        for iv in range(nv):
            for ix in range(nx):
                dtfion_tmp = MATHut.Deriv1(fion[:, ix, iv], dt,
                                           periodic=0, order=FD_order)
                dtfion[:, ix, iv] = dtfion_tmp[nb_ghost:-nb_ghost]

    dic_POPEop_VlasovLHS = {}
    dic_POPEop_VlasovLHS['timegrid'] = timegrid_V
    dic_POPEop_VlasovLHS['xgrid'] = xgrid
    dic_POPEop_VlasovLHS['vgrid'] = vgrid
    dic_POPEop_VlasovLHS['dxPhi'] = dxPhi
    dic_POPEop_VlasovLHS['v_dxfelec'] = v_dxfelec
    dic_POPEop_VlasovLHS['dvfelec'] = dvfelec
    dic_POPEop_VlasovLHS['dtfelec'] = dtfelec
    if fion_exist:
        dic_POPEop_VlasovLHS['v_dxfion'] = v_dxfion
        dic_POPEop_VlasovLHS['dvfion'] = dvfion
        dic_POPEop_VlasovLHS['dtfion'] = dtfion

    return dic_POPEop_VlasovLHS

# end def Compute_POPEop_forVlasovLHS


# --------------------------------------------------
#  Get (load if already exist or compute)
#   the operators appearing in Left Had Side of
#   Vlasov equations
# --------------------------------------------------
def Get_POPEop_forVlasovLHS(
        FD_order, # Finite Difference order
        POPE_read,
        POPE_h5filename):
    '''
    Load the LHS operators of the Vlasov equation saved in POPE_h5filename
    '''

    h5group_Vlasov = 'Vlasov/'
    if os.path.exists(POPE_h5filename):
        dic_POPEop_Vlasov_exist = H5ut.group_exist_in_hdf5(POPE_h5filename,
                                                           h5group_Vlasov)
        if dic_POPEop_Vlasov_exist:
            print('==> Loading of POPE operators for Vlasov in '+POPE_h5filename)
            dic_POPEop_Vlasov = H5ut.load_dict_from_hdf5(POPE_h5filename,
                                                         h5group_Vlasov)
            return dic_POPEop_Vlasov

    dic_POPEop_Vlasov = Compute_POPEop_forVlasovLHS(
        FD_order, # Finite Difference order
        POPE_read)

    print('==> Saving of POPE operators for Vlasov in '+POPE_h5filename)
    H5ut.save_dict_to_hdf5(dic_POPEop_Vlasov, h5group_Vlasov, POPE_h5filename)

    return dic_POPEop_Vlasov

# end def Get_POPEop_forVlasov_LHS


# --------------------------------------------------
#  Get (load if already exist or compute)
#   all the operators appearing in Vlasov equation
#   including RHS
# --------------------------------------------------
def Get_POPEop_forVlasov(
        FD_order, # Finite Difference order
        POPE_read,
        POPE_h5filename):
    '''
    Load operators of the Vlasov equation saved in POPE_h5filename
    '''

    dic_POPEop_VlasovLHS = Get_POPEop_forVlasovLHS(
        FD_order, # Finite Difference order
        POPE_read,
        POPE_h5filename)

    return [dic_POPEop_VlasovLHS]

# end def Get_POPEop_forVlasov
