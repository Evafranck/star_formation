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

from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import astropy 

# Define model
from astropy.modeling import Fittable1DModel, Parameter

from scipy.optimize import curve_fit

#plt.switch_backend('agg') 
#%matplotlib inline
#%config InlineBackend.figure_format='retina'


#####################################################################################################
titlelist = ['Low Resolution', 'Medium Resolution',  'High Resolution']
path = ['../low_master_iso', '../low_semenov_iso', '../low_evans_iso', '../low_federrath_iso']
path2 = ['../med_master_iso', '../med_semenov_iso', '../med_evans_iso', '../med_federrath_iso']
path3 = ['../high_master_iso', '../high_semenov_iso', '../high_evans_iso', '../high_federrath_iso']
simulation = ['Threshold-based model','Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']

pathlist = [path, path2, path3]
######################################################################################################

fig = plt.figure(figsize=(15,5))
gs0 = gridspec.GridSpec(1, 3, height_ratios = [1], width_ratios = [1, 1, 1])
gs0.update(hspace=0.00, wspace=0.00)

label = ['Threshold-based model', 'Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2012)']
color = ['blue', 'orange', 'green', 'red']
loc = [0.4, 0.6]
for n in range(3):
    ax = fig.add_subplot(gs0[n])
    ax.set_title(titlelist[n], fontsize = 13)
    ax.set_xlabel('time [Gyr]', fontsize = 12)
    if n==1:
        ax.set_yticklabels([])
    if n==0:
        ax.set_ylabel(r'Gas mass [log($M_{\rm sol}$)]', fontsize = 12)
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
            disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
            s_disk = s_all[disk]
            cold = f.LowPass('temp', '15000 K') # nur kaltes gas
            s = s_disk.g[cold]
            gas_mass.append(np.sum(s.g['mass']))
            time.append(s.properties['time'].in_units('Gyr'))
            
        # Sortiere die Indizes der x-Werte
        sorted_indices = sorted(range(len(time)), key=lambda i: time[i])

        # Verwende die sortierten Indizes, um x und y neu zu ordnen
        time_sorted = [time[i] for i in sorted_indices]
        gas_mass_sorted = [gas_mass[i] for i in sorted_indices]
        ax.plot(time_sorted,np.log10(gas_mass_sorted),color=color[k],lw=1, linestyle = '-', label = label[k])
        #ax.set_ylim(-0.1,4.45)
        #ax.set_xlim(0.1,1.6)

for n in range(2):
    if n==0: ax.legend(fontsize=10, loc='lower left',bbox_to_anchor=(-1.95, 0.05))

plt.savefig('gas_mass_time.pdf', bbox_inches='tight')
