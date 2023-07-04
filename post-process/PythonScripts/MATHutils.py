# SPDX-License-Identifier: MIT
'''
File containing mathematical tools
'''

import math
import numpy as np
from scipy import integrate

# *** COEFFICIENTS FOR FINITE DIFFERENCE ***
dic_coeff_FD = {}

# *** Finite differences for first and second derivative at 2th order ***
denom = 2.
dic_coeff_FD['[1,2],centered'] = np.array([-1., 0., 1.])/denom
dic_coeff_FD['[1,2],forward'] = np.array([-3., 4., -1.])/denom
dic_coeff_FD['[1,2],backward'] = np.array([1., -4., 3.])/denom

dic_coeff_FD['[2,2],centered'] = np.array([1., -2., 1.])
dic_coeff_FD['[2,2],forward'] = np.array([1., -2., 1.])
dic_coeff_FD['[2,2],backward'] = np.array([1., -2., 1.])

# *** Finite differences for first and second derivative at 4th order ***
denom = 12.
# - centered coefficients (-2,-1,0,1,2): Rk c_{i} = -c{-i} and c_0 = 0.
dic_coeff_FD['[1,4],centered'] = np.array([1., -8., 0., 8., -1.])/denom
# - forward and backward coefficients:
# 'forward' for (0,1,2,3,4) and 'backward' for (-4,-3,-2,-1,0)
dic_coeff_FD['[1,4],forward'] = np.array([-25., 48., -36., 16., -3.])/denom
dic_coeff_FD['[1,4],backward'] = np.array([3., -16., 36., -48., 25.])/denom
# - mixed coefficients 'mixed_p' for (-1,0,1,2,3) and 'mixed_m' for (-3,-2,-1,0,1)
dic_coeff_FD['[1,4],mixed_p'] = np.array([-3., -10., 18, -6., 1.])/denom
dic_coeff_FD['[1,4],mixed_m'] = np.array([-1., 6., -18, 10., 3.])/denom

denom = 12.
dic_coeff_FD['[2,4],centered'] = np.array([-1., 16., -30., 16., -1.])/denom
dic_coeff_FD['[2,4],forward'] = np.array([45., -154., 214., -156., 61., -10.])/denom
dic_coeff_FD['[2,4],backward'] = np.array([-10., 61., -156., 214., -154., 45.])/denom
dic_coeff_FD['[2,4],mixed_p'] = np.array([10., -15., -4., 14., -6., 1.])/denom
dic_coeff_FD['[2,4],mixed_m'] = np.array([1., -6., 14., -4., -15., 10.])/denom

# *** Finite differences for first and second derivative at 6th order ***
denom = 60.
dic_coeff_FD['[1,6],centered'] = np.array([-1., 9., -45., 0., 45., -9., 1.])/denom
dic_coeff_FD['[1,6],forward'] = np.array([-147., 360., -450.,
                                          400., -225., 72., -10.])/denom
dic_coeff_FD['[1,6],backward'] = np.array([10., -72., 225., -400.,
                                           450., -360., 147.])/denom
dic_coeff_FD['[1,6],mixed_p'] = np.array([-10., -77., 150.,
                                          -100., 50., -15., 2.])/denom
dic_coeff_FD['[1,6],mixed_m'] = np.array([-2., 15., -50.,
                                          100., -150, 77., 10.])/denom

denom = 180.
dic_coeff_FD['[2,6],centered'] = np.array([2., -27., 270.,
                                           -490., 270., -27., 2.])/denom
dic_coeff_FD['[2,6],forward'] = np.array([938., -4014., 7911., -9490.,
                                          7380., -3618., 1019., -126.])/denom
dic_coeff_FD['[2,6],backward'] = np.array([-126., 1019., -3618., 7380.,
                                           -9490., 7911., -4014., 938.])/denom
dic_coeff_FD['[2,6],mixed_p'] = np.array([126., -70., -486., 855.,
                                          -670., 324., -90., 11.])/denom
dic_coeff_FD['[2,6],mixed_m'] = np.array([11., -90., 324., -670.,
                                          855., -486., -70., 126.])/denom


# ---------------------------------------------------------------
#  nth derivative computed with finite difference
#    Input: F         = function to be derivate
#           dx        = step of the variable for derivative
#           periodic  = 1 if F is periodic
#                     = 0 otherwise (by default)
#           nth_deriv = nth derivative (1 or 2 (1 by default))
#           order     = error order (2,4 or 6 (4 by default))
#    Output: dnFdxn = nth derivative of F
# ---------------------------------------------------------------
def Deriv_nthorder(F, dx, periodic=0, nth_deriv=1, order=4):
    '''
    nth derivative computed with finite difference
       Input: F         = function to be derivate
              dx        = step of the variable for derivative
              periodic  = 1 if F is periodic
                        = 0 otherwise (by default)
              nth_deriv = nth derivative (1 or 2 (1 by default))
              order     = error order (2,4 or 6 (4 by default))
       Output: dnFdxn = nth derivative of F
    '''

    if nth_deriv not in [1,2]:
        raise ValueError(str(nth_deriv) +
                         "= unexpected derivative : should be 1 or 2")
    if order not in [2,4,6]:
        raise ValueError(str(order) + "= unexpected order : should be 2,4 or 6")

    nx = len(F)
    dnFdxn = np.zeros((nx))

    offset = order//2

    s_deriv_order = str([nth_deriv,order]).replace(" ", "")
    cc = dic_coeff_FD[s_deriv_order+',centered']
    cb = dic_coeff_FD[s_deriv_order+',backward']
    cf = dic_coeff_FD[s_deriv_order+',forward']
    if order > 2:
        cm_m = dic_coeff_FD[s_deriv_order+',mixed_m']
        cm_p = dic_coeff_FD[s_deriv_order+',mixed_p']

    # centered differences for dnFdxn[i] with offset<i<nx-offset
    #  --> example for 4th order :
    #   dnFdxn[i] = c[0](F[i-2]-F[i+2]) + c[1](F[i-1]-F[i+1])
    for i, coeff in enumerate(cc):
        dnFdxn[offset:nx-offset] = dnFdxn[offset:nx-offset] + coeff*F[i:nx-order+i]

    if periodic == 0:
        # forward differences for dnFdxn[0] and dnFdxn[nx-1] of 4th order
        #  dnFdxn[0] = cf_0*F[0] + cf_1*F[1] + cf_2*F[2] + cf_3*F[3] + cf_4*F[4]
        dnFdxn[0] = np.dot(cf, F[0:len(cf)])
        dnFdxn[nx-1] = np.dot(cb, F[nx-len(cb):])
        # mixed differences for dnFdxn[1] and dnFdxn[nx-2] of 4th order
        if order > 2:
            for indx in np.arange(offset-1):
                dnFdxn[indx+1] = np.dot(cm_p, F[indx:indx+len(cm_p)])
                dnFdxn[nx-2-indx] = np.dot(cm_m, F[nx-indx-len(cm_m):nx-indx])
    else:
        # boundary indexes to fill -> ex for 4th order: 0,1,nx-2,nx-1
        indx_bound = np.concatenate((np.arange(offset), nx-offset+np.arange(offset)))
        for indx in indx_bound:
            for i, coeff in enumerate(cc):
                dnFdxn[indx] = dnFdxn[indx] + coeff*F[np.fmod(indx+i-offset, nx)]

    dnFdxn = dnFdxn / math.pow(dx, nth_deriv)

    return dnFdxn

# end def Deriv_nthorder


# ---------------------------------------------------------------
#  First derivative computed with finite difference
# ---------------------------------------------------------------
def Deriv1(F, dx, periodic=0, order=4):
    '''
     First derivative computed with finite difference
       Input: F         = function to be derivate
              dx        = step of the variable for derivative
              periodic  = 1 if F is periodic
                        = 0 otherwise (by default)
              order     = error order (2,4 or 6 (4 by default))
       Output: dFdx = first derivative of F
    '''

    dFdx = Deriv_nthorder(F, dx, periodic=periodic, nth_deriv=1 , order=order)
    return dFdx

# end Deriv1


# ---------------------------------------------------------------
#  Second derivative computed with finite difference
# ---------------------------------------------------------------
def Deriv2(F, dx, periodic=0, order=4):
    '''
    Second derivative computed with finite difference
       Input: F         = function to be derivate
              dx        = step of the variable for derivative
              periodic  = 1 if F is periodic
                        = 0 otherwise (by default)
              order     = error order (2,4 or 6 (4 by default))
       Output: d2Fdx2 = second derivative of F
    '''

    d2Fdx2 = Deriv_nthorder(F, dx, periodic=periodic, nth_deriv=2, order=order)
    return d2Fdx2

# end Deriv2


# -------------------------------------------------
#  Personal FFT1D function
# -------------------------------------------------
def Fourier1D(F, x):
    ''' Personal FFT1D function '''

    x = np.asarray(x)
    assert x.ndim == 1
    dx = np.abs(x[1]-x[0])

    kx = (2.*np.pi/dx)*np.fft.fftfreq(x.shape[-1])
    TFF = np.fft.fft(F)

    kx = np.fft.fftshift(kx)
    TFF = np.fft.fftshift(TFF)

    return TFF, kx

# end def Fourier1D


# -------------------------------------------------
#  Personal inverse FFT1D function
# -------------------------------------------------
def iFourier1D(TFF, kx):
    ''' Personal inverse FFT1D function '''

    kx = np.asarray(kx)
    assert kx.ndim == 1

    nx = len(kx)
    dkx = np.abs(kx[1] - kx[0])
    x = np.linspace(0., 2.*np.pi/dkx, nx, endpoint=False)
    TFF = np.fft.ifftshift(TFF)
    F = np.real(np.fft.ifft(TFF))

    return F, x

# end def iFourier1D


# -------------------------------------------------
#  Personal FFT2D function
# -------------------------------------------------
def Fourier2D(F, x, y):
    ''' Personal FFT2D function '''

    x = np.asarray(x)
    assert x.ndim == 1
    dx = np.abs(x[1]-x[0])

    y = np.asarray(y)
    assert y.ndim == 1
    dy = np.abs(y[1]-y[0])

    kx = (2.*np.pi/dx)*np.fft.fftfreq(x.shape[-1])
    ky = (2.*np.pi/dy)*np.fft.fftfreq(y.shape[-1])
    TFF = np.fft.fft2(F)

    kx = np.fft.fftshift(kx)
    ky = np.fft.fftshift(ky)
    TFF = np.fft.fftshift(TFF)

    return TFF, kx, ky

# end def Fourier2D


# -------------------------------------------------
#  Personal inverse FFT2D function
# -------------------------------------------------
def iFourier2D(TFF, kx, ky):
    ''' Personal inverse FFT2D function '''

    kx = np.asarray(kx)
    assert kx.ndim == 1
    ky = np.asarray(ky)
    assert ky.ndim == 1

    nx = len(kx)
    dkx = np.abs(kx[1] - kx[0])
    x = np.linspace(0., 2.*np.pi/dkx, nx, endpoint=False)
    ny = len(ky)
    dky = np.abs(ky[1] - ky[0])
    y = np.linspace(0., 2.*np.pi/dky, ny, endpoint=False)

    TFF = np.fft.ifftshift(TFF)
    F = np.real(np.fft.ifft2(TFF))

    return F, x, y

# end def iFourier2D


# --------------------------------------------------
#  Computation of int F(v,x,t) dv
# --------------------------------------------------
def compute_intFdv(F, v):
    '''
    Computation of int F(v,x,t) dv
    '''

    [nt, nx, _] = np.shape(F)
    intFdv = np.zeros([nt, nx])
    for ix in np.arange(nx):
        for it in np.arange(nt):
            intFdv[it, ix] = integrate.simps(F[it, ix, :], v)

    return intFdv

# end def compute_intFdv


# --------------------------------------------------
#   Computation of PDF
# --------------------------------------------------
def Compute_PDF(data, theory_val):
    '''
    Computation of the PDF for the data set given in the variable 'data'
    using the theoretical value 'theory_val'
    '''

    # Computation of the error on the data set.
    # The commented line below does the same thing but
    # in linear scale

    error = np.log10(np.abs(data - theory_val))
    #error = data - theory_val

    # Computation of the foundamental parameters
    # of the error data set [min, max, mean, rms]

    min_err = np.min(error)
    max_err = np.max(error)
    mean = np.mean(error)
    rms = np.sqrt(np.mean((error-mean)*(error-mean)))

    # Binning operation on the error axis

    step = rms/100
    Np_p = int((max_err - mean)/step)+1  # Number of positive bins
    Np_m = int((mean-min_err)/step)+1    # Number of negative bins
    N = Np_m + Np_p                      # Number of total bins
    stat_min = mean - Np_m*step          # Minimum value of the error axis
    stat_max = mean + Np_p*step          # Maximum value of the error axis
    bins = np.arange(stat_min, stat_max, step) # Error axis

    # Building the NDF
    IND = (error - stat_min)/step
    IND = IND.astype(int) + 1
    NDF = np.zeros(N)

    for i in np.arange(len(IND)):
        if (IND[i] > 0 & IND[i] < N):
            a = IND[i]
            NDF[a-1] = NDF[a-1] + 1

    # Building the PDF
    PDF = NDF / np.sum(NDF)
    mean_PDF = np.mean(np.exp(bins*np.log(10))*PDF)
    #mean_PDF = np.sum(PDF*bins)/len(PDF)
    print('mean_PDF = '+str(mean_PDF))


    #VG#""" Plots"""
    #VG#mpp.plot(bins, NDF,'.')
    #VG#mpp.yscale('log')
    #VG#mpp.title('NDF',size=20)
    #VG#mpp.xlabel('be-1',size=20)
    #VG#mpp.savefig('NDF_be_dt2-5_log.png')
    #VG#mpp.show()
    #VG#
    #VG#mpp.plot(bins,PDF,'.')
    #VG#mpp.yscale('log')
    #VG#mpp.title('PDF', size = 20)
    #VG#mpp.xlabel('be-1',size=20)
    #VG#mpp.savefig('PDF_be_dt2-5_log.png')
    #VG#mpp.show()

# end def Compute_PDF


# ------------------------------------------------------
#  Find the lowest non-zero number in the array and
#  replace all zeroes by a number lower than the lowest
# ------------------------------------------------------
def replaceZeros(data):
    '''
    Find the lowest non-zero number in the array and
    replace all zeroes by a number lower than the lowest
    '''

    min_nonzero = np.min(data[np.nonzero(data)])
    data[data == 0] = min_nonzero*0.001

    return data

# end def replaceZeros
