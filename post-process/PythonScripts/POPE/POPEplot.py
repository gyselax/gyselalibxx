'''
Plot figures related to POPE method

  Author: V. Grandgirard
'''

# pylint: disable=invalid-name

# --- Matplotlib configuration ---
import matplotlib

import matplotlib.pyplot as plt
import numpy as np

import MATHutils as MATHut

matplotlib.rcParams['axes.formatter.limits'] = 3, -3
matplotlib.rcParams['figure.max_open_warning'] = 0

plt.viridis()
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

ticker_strformat = '%1.2f'

# --------------------------------------------------
# --> Used to plot POPE coefficients
# --------------------------------------------------
def my_plot_(fig, ax,
             x, y, f_xy,
             xlabel, ylabel, title):
    '''
    Definition of personal plot
    '''

    if len(y) == 1:
        ax.plot(x, np.reshape(f_xy, -1))
        ax.set_xlabel(xlabel)
    else:
        p = ax.pcolormesh(x, y, f_xy)
        fig.colorbar(p)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    ax.set_title(title)

# end def my_plot_


# --------------------------------------------------
# --> Plot the electrostatic potential and
#    the distribution functions
# --------------------------------------------------
def Plot_Phi_fdistribu(POPE_read,
                       time_acronym,
                       POPE_dir):
    '''
    Plot the electrostatic potential and
     the distribution functions
    '''

    xgrid = POPE_read.xgrid
    vgrid = POPE_read.vgrid
    timegrid = POPE_read.timegrid
    nb_ghosts = POPE_read.nb_ghosts
    Phi = POPE_read.Phi
    felec = POPE_read.felec
    if hasattr(POPE_read, 'fion'):
        fion = POPE_read.fion

    [nt, _, nv] = np.shape(felec)
    indx_v = int(nv/3)
    str_v = str(round(vgrid[indx_v], 4))

    fig = plt.figure(figsize=(18, 6))
    ax1 = fig.add_subplot(131)
    my_plot_(fig, ax1,
             xgrid, timegrid[nb_ghosts:nt-nb_ghosts],
             Phi[nb_ghosts:nt-nb_ghosts, :],
             'x', 't', 'Phi')

    ax2 = fig.add_subplot(132)
    my_plot_(fig, ax2,
             xgrid, timegrid[nb_ghosts:nt-nb_ghosts],
             felec[nb_ghosts:nt-nb_ghosts, :, indx_v],
             'x', 't', 'felec at v='+str_v)

    if hasattr(POPE_read, 'fion'):
        ax3 = fig.add_subplot(133)
        my_plot_(fig, ax3,
                 xgrid, timegrid[nb_ghosts:nt-nb_ghosts],
                 fion[nb_ghosts:nt-nb_ghosts, :, indx_v],
                 'x', 't', 'fion at v='+str_v)

    save_filename = POPE_dir+'/Phi_fdistribu_'+time_acronym+'.png'
    plt.savefig(save_filename)
    print('--> Saving of figure ' + save_filename)

# end def Plot_Phi_fdistribu


# --------------------------------------------------
# --> Plot the POPE coefficients for Poisson
#      equation
# --------------------------------------------------
def Plot_POPEcoeff_forPoisson(
        dic_POPEcoeff_Poisson,
        time_acronym,
        POPE_dir):
    '''
    Plot the POPE coefficients for Poisson equation
    '''

    xgrid = dic_POPEcoeff_Poisson['xgrid']
    timegrid = dic_POPEcoeff_Poisson['timegrid']
    qe = dic_POPEcoeff_Poisson['qe']
    qi = dic_POPEcoeff_Poisson['qi']

    [_, nx] = np.shape(qe)
    nx_bound = 2

    # ************ Plot POPE coeff for Poisson *****************
    fig = plt.figure(figsize=(12, 6))
    fig.suptitle("POPE coeff for Poisson solver")

    ax1 = fig.add_subplot(121)
    my_plot_(fig, ax1,
             xgrid[nx_bound:nx-nx_bound], timegrid, qe[:, nx_bound:nx-nx_bound],
             'x', 't', 'POPE qe')

    ax2 = fig.add_subplot(122)
    my_plot_(fig, ax2,
             xgrid[nx_bound:nx-nx_bound], timegrid, qi[:, nx_bound:nx-nx_bound],
             'x', 't', 'POPE qi')

    save_filename = POPE_dir+'/POPEcoeff_Poiss_'+time_acronym+'.png'
    plt.savefig(save_filename)
    print('--> Saving of figure ' + save_filename)

# end def Plot_POPEcoeff_forPoisson


# --------------------------------------------------
# --> Plot the POPE coefficients for Vlasov
#      equation
# --------------------------------------------------
def Plot_POPEcoeff_forVlasov(
        dic_POPEcoeff_Vlasov,
        time_acronym,
        POPE_dir):
    '''
    Plot the POPE coefficients for Vlasov equation
    '''

    timegrid = dic_POPEcoeff_Vlasov['timegrid']
    xgrid = dic_POPEcoeff_Vlasov['xgrid']
    vgrid = dic_POPEcoeff_Vlasov['vgrid']
    ae = dic_POPEcoeff_Vlasov['ae']
    be = dic_POPEcoeff_Vlasov['be']
    if 'ce' in dic_POPEcoeff_Vlasov.keys():
        ce = dic_POPEcoeff_Vlasov['ce']
        de = dic_POPEcoeff_Vlasov['de']
    if 'ai' in dic_POPEcoeff_Vlasov.keys():
        ai = dic_POPEcoeff_Vlasov['ai']
        bi = dic_POPEcoeff_Vlasov['bi']
    if 'ci' in dic_POPEcoeff_Vlasov.keys():
        ci = dic_POPEcoeff_Vlasov['ci']
        di = dic_POPEcoeff_Vlasov['di']

    [_, nx, nv] = np.shape(ae)
    nx_bound = 2
    indx_v = int(nv/3)
    str_v = str(round(vgrid[indx_v], 4))

    # ************ Plot POPE coeff for Vlasov for electrons *****************
    fig = plt.figure(figsize=(12, 12))
    fig.suptitle("POPE coeff for Vlasov solver for electrons")

    ax1 = fig.add_subplot(221)
    ax1_title = 'POPE ae at v = '+str_v
    my_plot_(fig, ax1,
             xgrid[nx_bound:nx-nx_bound], timegrid,
             ae[:, nx_bound:nx-nx_bound, indx_v],
             'x', 't', ax1_title)

    ax2 = fig.add_subplot(222)
    ax2_title = 'POPE be at v = '+str_v
    my_plot_(fig, ax2,
             xgrid[nx_bound:nx-nx_bound], timegrid,
             be[:, nx_bound:nx-nx_bound, indx_v],
             'x', 't', ax2_title)

    if 'ce' in dic_POPEcoeff_Vlasov.keys():
        ax3 = fig.add_subplot(223)
        ax3_title = 'POPE ce at v = '+str_v
        my_plot_(fig, ax3,
                 xgrid[nx_bound:nx-nx_bound], timegrid,
                 ce[:, nx_bound:nx-nx_bound, indx_v],
                 'x', 't', ax3_title)

    if 'de' in dic_POPEcoeff_Vlasov.keys():
        ax4 = fig.add_subplot(224)
        ax4_title = 'POPE de at v = '+str_v
        my_plot_(fig, ax4,
                 xgrid[nx_bound:nx-nx_bound], timegrid,
                 de[:, nx_bound:nx-nx_bound, indx_v],
                 'x', 't', ax4_title)

    save_filename = POPE_dir+'/POPEcoeff_elec_'+time_acronym+'.png'
    plt.savefig(save_filename)
    print('--> Saving of figure '+save_filename)

    # ************ Plot POPE coeff for Vlasov for ions *****************
    if 'ai' in dic_POPEcoeff_Vlasov.keys():
        fig = plt.figure(figsize=(12, 12))
        fig.suptitle("POPE coeff for Vlasov solver for ions")

        ax1 = fig.add_subplot(221)
        ax1_title = 'POPE ai at v = '+str_v
        my_plot_(fig, ax1,
                 xgrid[nx_bound:nx-nx_bound], timegrid,
                 ai[:, nx_bound:nx-nx_bound, indx_v],
                 'x', 't', ax1_title)

        ax2 = fig.add_subplot(222)
        ax2_title = 'POPE bi at v = '+str_v
        my_plot_(fig, ax2,
                 xgrid[nx_bound:nx-nx_bound], timegrid,
                 bi[:, nx_bound:nx-nx_bound, indx_v],
                 'x', 't', ax2_title)

        if 'ci' in dic_POPEcoeff_Vlasov.keys():
            ax3 = fig.add_subplot(223)
            ax3_title = 'POPE ci at v = '+str_v
            my_plot_(fig, ax3,
                     xgrid[nx_bound:nx-nx_bound], timegrid,
                     ci[:, nx_bound:nx-nx_bound, indx_v],
                     'x', 't', ax3_title)

        if 'di' in dic_POPEcoeff_Vlasov.keys():
            ax4 = fig.add_subplot(224)
            ax4_title = 'POPE di at v = '+str_v
            my_plot_(fig, ax4,
                     xgrid[nx_bound:nx-nx_bound], timegrid,
                     di[:, nx_bound:nx-nx_bound, indx_v],
                     'x', 't', ax4_title)

        save_filename = POPE_dir+'/POPEcoeff_ion_'+time_acronym+'.png'
        plt.savefig(save_filename)
        print('--> Saving of figure '+save_filename)

# end def Plot_POPEcoeff_forVlasov


# --------------------------------------------------
# --> Plot the iPOPE coefficients for Poisson
# --------------------------------------------------
def Plot_iPOPEcoeff_forPoisson(
        dic_iPOPEcoeff_Poisson,
        time_acronym,
        POPE_dir):
    '''
    Plot the iPOPE coefficients for Poisson
    '''

    xgrid = dic_iPOPEcoeff_Poisson['xgrid']
    iPOPE_index_q = dic_iPOPEcoeff_Poisson['iPOPE_index_q']
    iPOPE_error_q = dic_iPOPEcoeff_Poisson['iPOPE_error_q']
    iPOPE_residu_q = dic_iPOPEcoeff_Poisson['iPOPE_residu_q']

    nx = len(iPOPE_index_q)
    nx_bound = 2

    fig = plt.figure(figsize=(18, 6))
    fig.suptitle("iPOPE coeff for Poisson solver")

    ax1 = fig.add_subplot(131)
    ax1.plot(xgrid[nx_bound:nx-nx_bound],
             iPOPE_index_q[nx_bound:nx-nx_bound])
    ax1.set_xlabel('x')
    ax1.set_title('iPOPE index=-log10|<Coeff_q>|')

    ax2 = fig.add_subplot(132)
    abs_iPOPE_error_q = MATHut.replaceZeros(np.abs(iPOPE_error_q))
    ax2.plot(xgrid[nx_bound:nx-nx_bound],
             np.log10(abs_iPOPE_error_q[nx_bound:nx-nx_bound]))
    ax2.set_xlabel('x')
    ax2.set_title('log10|<error_q>|')

    ax3 = fig.add_subplot(133)
    ax3.semilogy(xgrid[nx_bound:nx-nx_bound],
                 np.abs(iPOPE_residu_q[nx_bound:nx-nx_bound]))
    ax3.set_xlabel('x')
    ax3.set_title('<residu_q>')

    save_filename = POPE_dir+'/iPOPEcoeff_Poiss_'+time_acronym+'.png'
    plt.savefig(save_filename)
    print('--> Saving of figure ' + save_filename)

# end def Plot_iPOPEcoeff_forPoisson


# --------------------------------------------------
# --> Plot the iPOPE coefficients for
#      Vlasov equation
# --------------------------------------------------
def Plot_iPOPEcoeff_forVlasov(
        dic_iPOPEcoeff_Vlasov,
        time_acronym,
        POPE_dir):
    '''
    Plot the iPOPE coefficients for Vlasov equation
    '''

    xgrid = dic_iPOPEcoeff_Vlasov['xgrid']
    vgrid = dic_iPOPEcoeff_Vlasov['vgrid']
    iPOPE_index_e = dic_iPOPEcoeff_Vlasov['iPOPE_index_e']
    iPOPE_error_e = dic_iPOPEcoeff_Vlasov['iPOPE_error_e']
    iPOPE_residu_e = dic_iPOPEcoeff_Vlasov['iPOPE_residu_e']
    if 'iPOPE_index_i' in list(dic_iPOPEcoeff_Vlasov):
        fion_exist = True
        iPOPE_index_i = dic_iPOPEcoeff_Vlasov['iPOPE_index_i']
        iPOPE_error_i = dic_iPOPEcoeff_Vlasov['iPOPE_error_i']
        iPOPE_residu_i = dic_iPOPEcoeff_Vlasov['iPOPE_residu_i']
    else:
        fion_exist = False

    [nx, nv] = np.shape(iPOPE_index_e)
    nx_bound = 2
    nv_bound = 8

    # ************ Plot iPOPE analysis for electrons *****************
    fig = plt.figure(figsize=(18, 6))
    fig.suptitle("iPOPE coeff for Vlasov for electrons")

    ax1 = fig.add_subplot(131)
    p1 = ax1.pcolormesh(xgrid[nx_bound:nx-nx_bound],
                        vgrid[nv_bound:nv-nv_bound],
                        iPOPE_index_e[nx_bound:nx-nx_bound,
                                      nv_bound:nv-nv_bound].T)
    fig.colorbar(p1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('v')
    ax1.set_title('iPOPE index=-log10|<Coeff_e>|')

    ax2 = fig.add_subplot(132)
    abs_iPOPE_error_e = MATHut.replaceZeros(np.abs(iPOPE_error_e))
    p2 = ax2.pcolormesh(xgrid[nx_bound:nx-nx_bound],
                        vgrid[nv_bound:nv-nv_bound],
                        np.log10(abs_iPOPE_error_e[nx_bound:nx-nx_bound,
                                                   nv_bound:nv-nv_bound]).T)
    fig.colorbar(p2)
    ax2.set_xlabel('x')
    ax2.set_ylabel('v')
    ax2.set_title('log10|<error_e>|')

    ax3 = fig.add_subplot(133)
    p3 = ax3.pcolormesh(xgrid[nx_bound:nx-nx_bound],
                        vgrid[nv_bound:nv-nv_bound],
                        iPOPE_residu_e[nx_bound:nx-nx_bound,
                                       nv_bound:nv-nv_bound].T)
    fig.colorbar(p3)
    ax3.set_xlabel('x')
    ax3.set_ylabel('v')
    ax3.set_title('<residu_e>')

    save_filename = POPE_dir+'/iPOPEcoeff_elec_'+time_acronym+'.png'
    plt.savefig(save_filename)
    print('--> Saving of figure '+save_filename)

    # ************ Plot iPOPE analysis for ions *****************
    if fion_exist:
        fig = plt.figure(figsize=(18, 6))
        fig.suptitle("iPOPE coeff for ions")

        ax1 = fig.add_subplot(131)
        p1 = ax1.pcolormesh(xgrid[nx_bound:nx-nx_bound],
                            vgrid[nv_bound:nv-nv_bound],
                            iPOPE_index_i[nx_bound:nx-nx_bound,
                                          nv_bound:nv-nv_bound].T)
        fig.colorbar(p1)
        ax1.set_xlabel('x')
        ax1.set_ylabel('v')
        ax1.set_title('iPOPE index=-log10|<Coeff_i>|')

        ax2 = fig.add_subplot(132)
        abs_iPOPE_error_i = MATHut.replaceZeros(np.abs(iPOPE_error_i))
        p2 = ax2.pcolormesh(xgrid[nx_bound:nx-nx_bound],
                            vgrid[nv_bound:nv-nv_bound],
                            np.log10(abs_iPOPE_error_i[nx_bound:nx-nx_bound,
                                                       nv_bound:nv-nv_bound].T))
        fig.colorbar(p2)
        ax2.set_xlabel('x')
        ax2.set_ylabel('v')
        ax2.set_title('log10|<error_i>|')

        ax3 = fig.add_subplot(133)
        p3 = ax3.pcolormesh(xgrid[nx_bound:nx-nx_bound],
                            vgrid[nv_bound:nv-nv_bound],
                            iPOPE_residu_i[nx_bound:nx-nx_bound,
                                           nv_bound:nv-nv_bound].T)
        fig.colorbar(p3)
        ax3.set_xlabel('x')
        ax3.set_ylabel('v')
        ax3.set_title('<residu_i>')

        save_filename = POPE_dir+'/iPOPEcoeff_ion_'+time_acronym+'.png'
        plt.savefig(save_filename)
        print('--> Saving of figure '+save_filename)

# end def Plot_iPOPEcoeff_forVlasov
