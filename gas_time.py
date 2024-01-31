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

model = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan']

fig = plt.figure(figsize=(10,10))

label = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r' $\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)', 'Platzhalter']
color = ['blue', 'orange', 'green', 'red', 'purple', 'black']
plt.title('Gas mass over time', fontsize = 13)
plt.xlabel('time [Gyr]', fontsize = 12)
plt.ylabel(r'Gas mass [log($M_{\rm sol}$)]', fontsize = 12)

for k, sim in enumerate(model):
    
    gas_mass = []
    n = []
    time = []

    simname = glob.glob('../'+sim+'/halo.0????')
    print(simname[0])
    #if not os.path.isfile(sim+'_surf_den_test.dat'):
    for name in simname:
        s_all = pb.load(name)
        pb.analysis.angmom.faceon(s_all)
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
    plt.plot(time_sorted,np.log10(gas_mass_sorted),color=color[k],lw=1, linestyle = '-', label = label[k])
    #plt.ylim(-0.1,4.45)
    #plt.xlim(0.1,1)

plt.legend(fontsize=10)#, loc='lower left',bbox_to_anchor=(-1.95, 0.05))

plt.savefig('gas_mass_time.pdf', bbox_inches='tight')
