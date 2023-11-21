# function definitions to fit surface density

import numpy as np
import matplotlib.pylab as plt
import pynbody as pb
import os, pickle
import glob
from pathlib import Path
import matplotlib.gridspec as gridspec
import scipy
import pynbody.filt as f



#####################################################################################################
titlelist = ['Gas mass over time'+'\n'+ '(Low Resolution)', 'Gas mass over time'+'\n'+ '(Medium Resolution)']
path = ['../low_master_iso', '../low_semenov_iso', '../low_evans_iso', '../low_federrath_iso']
path2 = ['../med_master_iso', '../med_semenov_iso', '../med_evans_iso', '../med_federrath_iso']
path3 = ['../high_master_iso', '../high_semenov_iso', '../high_evans_iso', '../high_federrath_iso']
simulation = ['Threshold-based model','Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']

pathlist = [path, path2, path3]
######################################################################################################

fig = plt.figure(figsize=(10,10))
gs0 = gridspec.GridSpec(1, 1)
gs0.update(hspace=0.00, wspace=0.00)


label = ['Threshold-based model', 'Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']
color = ['blue', 'orange', 'green', 'red']
loc = [0.4, 0.6]
for n in range(1):
    ax = fig.add_subplot(gs0[n])
    #ax.text(0.7,loc[n],r'Gas mass',fontsize=12,color='black',horizontalalignment='center',verticalalignment='center')
    ax.set_title(titlelist[n], fontsize = 16)
    ax.set_xlabel('time [Gyr]', fontsize = 14)
    if n==1:
        ax.set_yticklabels([])
    if n==0:
        ax.set_ylabel('Gas mass [kpc]', fontsize = 14)
    for k, sim in enumerate(pathlist[n]):
        
        gas_mass = []
        n = []
        time = []
 
        simname = glob.glob(sim+'/'+'*.0????')
        print(simname[0])
        #if not os.path.isfile(sim+'_surf_den_test.dat'):
        for name in simname:
            s_all = pb.load(sim + '/' + name)
            s_all.physical_units()
            print(s_all.g['mass'].sum())
            disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
            s_disk = s_all[disk]
            cold = f.LowPass('temp', '15000 K') # nur kaltes gas
            s = s_disk.g[cold]
            gas_mass.append(s['mass'].sum())
            time.append(s.properties['time'].in_units('Gyr'))
        
        print(time)
        print(gas_mass)
        ax.plot(time,gas_mass,color=color[k], linestyle = '-', label = label[k])
        #ax.set_ylim(-0.1,4.45)
        #ax.set_xlim(0.1,1.6)

plt.legend(fontsize=12)

plt.savefig('gas_mass_time_both.pdf', bbox_inches='tight')
