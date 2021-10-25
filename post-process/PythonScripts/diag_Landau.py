# Using :
#  python ../post-process/PythonScripts/diag_Landau.py
#  python ../post-process/PythonScripts/diag_Landau.py --dir . --save
#  python ../post-process/PythonScripts/diag_Landau.py --dir . --no-save

from argparse import ArgumentParser
import glob
import matplotlib.pyplot as plt
import numpy as np
import sys

import HDF5utils as H5ut
import MATHutils as MATHut
import READutils as READut


def plot_Landau_damping(Phi, timegrid, time_end=None, dirname=".", plotfig=False, validate=False):
    """
    Plot the time evolution of abs(Phi) to compute the Landau damping rate and frequency
    """

    [nt, nx] = np.shape(Phi)

    irp = int(nx / 2)
    Phi_rp = Phi[:, irp]
    absPhi_rp = np.fabs(Phi_rp)
    if (time_end is None):
        nt_end = nt

    # Computation of the damping rate
    serie_max_indx = [0]
    serie_max_t = [timegrid[0]]
    serie_max_phi = [absPhi_rp[0]]
    for it in np.r_[0:nt_end - 1]:
        if ((absPhi_rp[it] > absPhi_rp[it - 1]) & (absPhi_rp[it] > absPhi_rp[it + 1])):
            serie_max_indx.append(it)
            serie_max_t.append(timegrid[it])
            serie_max_phi.append(absPhi_rp[it])
    indx_fit = 1
    icorr = 0
    for i in np.r_[1:len(serie_max_t) - 3]:
        correlation = np.corrcoef(serie_max_t[indx_fit:len(serie_max_t) - i],
                                  np.log(serie_max_phi[indx_fit:len(serie_max_t) - i]))[0, 1]
        correlation = correlation * correlation
        if (correlation > 0.999):
            icorr = i
        break
    poly_resu = np.polyfit(serie_max_t[indx_fit:len(serie_max_t) - icorr],
                           np.log(serie_max_phi[indx_fit:len(serie_max_t) - icorr]),
                           1, full='True')
    [p_fit, s_fit] = poly_resu[0]
    damping_rate = p_fit
    damping_rate_theory = -0.153
    damping_rate_relerror = abs(damping_rate - damping_rate_theory) / abs(damping_rate_theory)
    print("damping rate = ", damping_rate)
    print("--> rel. error  = ", damping_rate_relerror)
    if validate and damping_rate_relerror > 0.05:
        sys.exit(1)

    if (plotfig):
        plt.figure()
        plt.semilogy(timegrid[0:nt_end + 1], absPhi_rp[0:nt_end + 1])
        plt.semilogy(serie_max_t[:], serie_max_phi[:], 'o', markersize=8)
        plt.semilogy(timegrid[0:nt_end],
                     np.exp(s_fit) * np.exp(p_fit * timegrid[0:nt_end]),
                     'r', lw=3)
        plt.ylabel('abs(Phi(middle_box))', fontsize=16)
        plt.xlabel('time', fontsize=16)
        plt.title('damping_rate = ' + str(round(damping_rate, 4)) +
                  " ; theory = " + str(damping_rate_theory), fontsize=16)

        time_acronym = '_t' + str(round(timegrid[0], 4)
                                  ) + 'to' + str(round(timegrid[nt_end - 1], 4))
        savefig_filename = dirname + '/dampingrate' + time_acronym + '.png'
        print("Save figure in " + savefig_filename)
        plt.savefig(savefig_filename)

    # Computation of the frequency
    [FFTPhi_tt, freq] = MATHut.Fourier1D(Phi_rp, timegrid)
    abs_FFTPhi_tt = np.abs(FFTPhi_tt)
    frequency = np.abs(freq[np.max(np.nonzero(
        abs_FFTPhi_tt == np.max(abs_FFTPhi_tt)))])
    frequency_theory = 1.4156
    frequency_relerror = abs(frequency - frequency_theory) / frequency_theory
    print("frequency   = ", frequency)
    print(" --> rel. error = ", frequency_relerror)
    if validate and frequency_relerror > 0.05:
        sys.exit(1)

    if (plotfig):
        plt.figure()
        plt.plot(freq, FFTPhi_tt)
        plt.plot([frequency, frequency], [min(abs(FFTPhi_tt)), max(abs(FFTPhi_tt))],
                 'r-', lw=3, label='Code')
        plt.plot([frequency_theory, frequency_theory],
                 [min(abs(FFTPhi_tt)), max(abs(FFTPhi_tt))],
                 'g--', lw=3, label='Theory=' + str(frequency_theory))
        plt.ylabel('FFT(Phi(t))', fontsize=16)
        plt.xlabel('Freq', fontsize=16)
        plt.legend()
        plt.title('omega = ' + str(round(frequency, 4)), fontsize=16)
        time_acronym = '_t' + str(round(timegrid[0], 4)
                                  ) + 'to' + str(round(timegrid[nt_end - 1], 4))
        savefig_filename = dirname + '/frequency' + time_acronym + '.png'
        print("Save figure in " + savefig_filename)
        plt.savefig(savefig_filename)

    return [damping_rate, damping_rate_relerror, frequency, frequency_relerror]


# --------------------------------------------------
# MAIN
# --------------------------------------------------

def main():
    parser = ArgumentParser(description="Analysis of VOICEXX results")
    parser.add_argument('--dir', type=str, help="Folder containing the hdf5 files")
    save_group = parser.add_mutually_exclusive_group(required=False)
    save_group.add_argument('--save', dest='save', action='store_true', help="Save figures")
    save_group.add_argument(
        '--no-save',
        dest='save',
        action='store_false',
        help="Don't save figures")
    parser.set_defaults(save=None)
    args = parser.parse_args()

    if args.dir:
        VOICEXX_dir = args.dir
        if not READut.valid_directory(VOICEXX_dir):
            print('==> There are no VOICEXX results in this directory')
            sys.exit(1)
        else:
            print("Reading from directory : ", VOICEXX_dir)
    else:
        VOICEXX_dir = READut.Ask_directory()

    vxx_res = READut.Read_VOICEXX_results(VOICEXX_dir)
    gridx = vxx_res.MeshX
    dx = gridx[1] - gridx[0]

    [damping_rate, damping_rate_relerror, frequency, frequency_relerror] = \
        plot_Landau_damping(vxx_res.electric_potential,
                            vxx_res.time_saved,
                            dirname=VOICEXX_dir,
                            plotfig=args.save)


if __name__ == "__main__":
    main()
