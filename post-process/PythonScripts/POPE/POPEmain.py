'''
file : POPEmain.py
date : 20/07/2016

Main program for VOICE POPE verification
'''

# pylint: disable=invalid-name
# pylint: disable=import-error

import os
import pathlib

import numpy as np

import READutils as READut
import POPEread as POPEr
import POPEoperators as POPEop
import POPEdiag as POPEd
import POPEdiag_iPOPE as iPOPEd
import POPEplot as POPEpl
import PYutils as PYut

# --- Test if the PYTHONPATH has been correctly defined ---

POPE_python = pathlib.Path(__file__).parent.absolute()
os.sys.path.insert(1, POPE_python)

def POPE_main(VOICE_resu,
              FD_order,    # Finite Difference order
              ifirst_diag,
              nb_diag):
    '''
    POPE Analysis at a certain time
    '''

    # ---> number of existing diagnostics
    nb_files = len(VOICE_resu.time_saved)

    nb_ghosts = int(FD_order/2)+1

    ifirst_diag = ifirst_diag - nb_ghosts
    ifirst_diag = max(nb_ghosts, ifirst_diag)
    ifirst_diag = min(ifirst_diag, nb_files-nb_diag-2*nb_ghosts)
    iend_diag = ifirst_diag + nb_diag + 2*nb_ghosts

    # ---> Read VOICE results
    POPE_read = POPEr.Read_POPE_results(
        VOICE_resu, ifirst_diag, iend_diag)
    POPE_read.nb_ghosts = nb_ghosts

    POPE_dir = VOICEXX_resu.dirname
    time_acronym = 't'+str(round(POPE_read.timegrid[nb_ghosts], 4)) + \
        '_nbdiag'+str(nb_diag)
    POPE_dir = POPE_dir+'/POPE_analysis'
    if not os.path.exists(POPE_dir):
        print('Create the folder ', POPE_dir)
        os.system('mkdir '+POPE_dir)
    POPE_h5filename = POPE_dir+'/POPE_res_'+time_acronym+'.h5'

    # ---> Plot Phi and the distribution functions
    POPEpl.Plot_Phi_fdistribu(POPE_read, time_acronym, POPE_dir)

    # -------------------------------------------------
    #  Compute POPE and iPOPE analysis
    #   for Poisson equation
    # -------------------------------------------------
    # ---> Compute POPE operator for Poisson
    dic_POPEop_Poisson = POPEop.Get_POPEop_forPoisson(
        FD_order, POPE_read, POPE_h5filename)
    # ---> Compute and plot POPE coefficients for Poisson
    dic_POPEcoeff_Poisson = POPEd.Get_POPEcoeff_forPoisson(
        POPE_read, dic_POPEop_Poisson, POPE_h5filename)
    POPEd.ComputeGaussian_POPE_coeff(dic_POPEcoeff_Poisson)
    POPEpl.Plot_POPEcoeff_forPoisson(
        dic_POPEcoeff_Poisson, time_acronym=time_acronym, POPE_dir=POPE_dir)
    # ---> Compute and plot iPOPE coefficients for Poisson
    dic_iPOPEcoeff_Poisson = iPOPEd.Get_iPOPEcoeff_forPoisson(
        POPE_read, dic_POPEop_Poisson, POPE_h5filename)
    POPEpl.Plot_iPOPEcoeff_forPoisson(
        dic_iPOPEcoeff_Poisson, time_acronym=time_acronym, POPE_dir=POPE_dir)

    # -------------------------------------------------
    #  Compute POPE and iPOPE analysis
    #   for Vlasov equation
    # -------------------------------------------------
    # ---> Compute POPE operator for Vlasov
    [dic_POPEop_Vlasov] = POPEop.Get_POPEop_forVlasov(
        FD_order, POPE_read, POPE_h5filename)
    # ---> Compute and plot POPE coefficients for Vlasov
    dic_POPEcoeff_Vlasov = POPEd.Get_POPEcoeff_forVlasov(
        POPE_read, dic_POPEop_Vlasov, POPE_h5filename)
    POPEd.ComputeGaussian_POPE_coeff(dic_POPEcoeff_Vlasov)
    POPEpl.Plot_POPEcoeff_forVlasov(
        dic_POPEcoeff_Vlasov,
        time_acronym=time_acronym,
        POPE_dir=POPE_dir)
    # ---> Compute and plot iPOPE coefficients for Vlasov
    dic_iPOPEcoeff_Vlasov = iPOPEd.Get_iPOPEcoeff_forVlasov(
        POPE_read, dic_POPEop_Vlasov, POPE_h5filename)
    POPEpl.Plot_iPOPEcoeff_forVlasov(
        dic_iPOPEcoeff_Vlasov, time_acronym=time_acronym, POPE_dir=POPE_dir)

    # ---> Save all the .png figures in one .pdf file
    POPE_PDFfilename = POPE_h5filename.replace('.h5', '.pdf')
    print('--> Saving of ' + POPE_PDFfilename)
    print(' ')
    convert_cmd = 'convert '+POPE_dir+'/*'+time_acronym+'.png ' + POPE_PDFfilename
    os.system(convert_cmd)


# -----------------------------------------------------------
#   MAIN PROGRAM
# -----------------------------------------------------------
if __name__ == "__main__":

    VOICEXX_dir = READut.Ask_directory()
    VOICEXX_resu = READut.Read_VOICEXX_results(VOICEXX_dir)

    FD_order_ = 6
    nb_diag_ = 10
    nb_time_ = 3
    time_saved = VOICEXX_resu.time_saved
    time_list = np.linspace(time_saved[0], time_saved[-1], nb_time_)
    for time_diag in time_list:
        ifirst_diag_ = PYut.search_pos(time_diag, time_saved)
        POPE_main(VOICEXX_resu, FD_order_, ifirst_diag_, nb_diag_)
