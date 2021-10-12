import glob;
import matplotlib.pyplot as plt
import numpy as np

import HDF5utils as H5ut

#--------------------------------------------------
# Plot the time evolution of log(abs(Phi))
#--------------------------------------------------
def plot_Philogabs(Phi,timegrid,
                   time_end=None,dirname=''):

    [nt,nx]=np.shape(Phi)
    absPhi = np.fabs(Phi)
    logPhi = np.log(absPhi)
    if (time_end==None):
        nt_end = nt

    plt.figure()
    plt.plot(timegrid[0:nt_end+1],logPhi[0:nt_end+1,int(nx/2)])
    plt.ylabel('log(abs(Phi))',fontsize=16)
    plt.xlabel('time',fontsize=16)
    plt.title('Potential middle box',fontsize=18)

    time_acronym='_t'+str(round(timegrid[0],4))+'to'+str(round(timegrid[nt_end-1],4))
    savefig_filename = 'logPhi'+time_acronym+'.png'
    plt.savefig(savefig_filename)
    plt.show()


#--------------------------------------------------
# MAIN
#--------------------------------------------------
ResuFileList = sorted(glob.glob('VOICEXX_[0-9]*.h5'))
nbFiles = np.size(ResuFileList)

Hinit = H5ut.loadHDF5('VOICEXX_initstate.h5');
gridx = Hinit.MeshX
dx = gridx[1]-gridx[0]

H = H5ut.loadHDF5(ResuFileList)

Phi_tx = - np.gradient(H.efield,dx,axis=1)
plot_Philogabs(Phi_tx,H.time_saved)

