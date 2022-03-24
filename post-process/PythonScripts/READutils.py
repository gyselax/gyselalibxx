# SPDX-License-Identifier: MIT

#--------------------------------------------------
#  file : READutils.py
#  date : 2016/04/12
#
#  Reading of the VOICEXX files
#--------------------------------------------------

""" Tools for accessing results folders
"""
import glob
import os

import HDF5utils as H5ut

#===============================================
def valid_directory(ResuDir):
    """ Check if result directory is valid
    """
    cache_file  = os.path.join(ResuDir, 'combined_diag.h5')
    find_fnames = os.path.join(ResuDir,'VOICEXX_[0-9]*.h5')
    find_fnames = glob.glob(find_fnames)
    nb_fnames   = len(find_fnames)
    return nb_fnames!=0 or os.path.exists(cache_file)

#===============================================
def Ask_directory():
    """ Ask for the result directory
    """
    ResuDir = ''
    while True:
        ResuDir = input("Result directory [default = current]: ? ")
        if ( ResuDir==''):
            ResuDir = str(os.getcwd())
        #end if
        if valid_directory(ResuDir):
            break
        else:
            print('==> There are no VOICEXX results in this directory')

    return ResuDir

#end def
#-----------------------------------------------


#===============================================
def Read_VOICEXX_results(ResuDir, cache = True):
    """ Read the result directory
    """
    cache_file = os.path.join(ResuDir, 'combined_diag.h5')
    ResuFileList = sorted(glob.glob(ResuDir+'/VOICEXX_[0-9]*.h5'))

    if cache and os.path.exists(cache_file) and \
            (len(ResuFileList)==0 or \
             os.path.getmtime(cache_file) > os.path.getmtime(ResuFileList[-1])):
        H = H5ut.loadHDF5(cache_file)
        H.dirname = ResuDir
        return H

    H = H5ut.loadHDF5(ResuFileList)
    Hinit = H5ut.loadHDF5(ResuDir+'/VOICEXX_initstate.h5')
    H.append(Hinit)
    H.dirname = ResuDir

    if cache:
        H5ut.write_results(cache_file, H)

    return H

#end def Read_VOICE_results

