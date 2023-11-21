# function definitions to fit surface density

import numpy as np
import matplotlib.pylab as plt
import pynbody as pb
import os, pickle
import glob
from pathlib import Path
import matplotlib.gridspec as gridspec
import scipy
import pynbody.filt as filt

from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import astropy 

# Define model
from astropy.modeling import Fittable1DModel, Parameter

from scipy.optimize import curve_fit

#plt.switch_backend('agg') 
#%matplotlib inline
#%config InlineBackend.figure_format='retina'

import seaborn as sns
plt.rcParams['figure.figsize'] = (15, 10)

sns.set_style('ticks')
#sns.set_style('darkgrid')
sns.set_context("talk",font_scale=1,rc={"lines.linewidth": 2,"axes.linewidth": 2})

plt.rc('axes', linewidth=1)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 3
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['axes.edgecolor'] = 'k'#'gray'
#plt.rcParams['axes.grid'] = True
#plt.rcParams['grid.color'] = 'lightgray'
#plt.rcParams['grid.linestyle'] = 'dashed' #dashes=(5, 1)
plt.rcParams['lines.dashed_pattern'] = 10, 3
plt.rcParams['grid.linewidth'] = 1.5
#plt.rcParams['axes.facecolor'] = 'whitesmoke'
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fancybox'] = True
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.shadow'] = False
plt.rcParams['legend.edgecolor'] = 'darkgray'
plt.rcParams['patch.linewidth'] = 1


#####################################################################################################
titlelist = ['Gas mass over time'+'\n'+ '(Low Resolution)', 'Gas mass over time'+'\n'+ '(Medium Resolution)']
path = ['../low_master_iso', '../low_semenov_iso', '../low_evans_iso', '../low_federrath_iso']
path2 = ['../med_master_iso', '../med_semenov_iso', '../med_evans_iso', '../med_federrath_iso']
path3 = ['../high_master_iso', '../high_semenov_iso', '../high_evans_iso', '../high_federrath_iso']
simulation = ['Threshold-based model','Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']

pathlist = [path, path2, path3]
######################################################################################################

fig = plt.figure(figsize=(10,5))
gs0 = gridspec.GridSpec(1, 2)
gs0.update(hspace=0.00, wspace=0.00)


label = ['Threshold-based model', 'Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']
color = ['blue', 'orange', 'green', 'red']
loc = [0.4, 0.6]
for n in range(3):
    ax = fig.add_subplot(gs0[n])
    ax.text(0.7,loc[n],r'Gas mass',fontsize=12,color='black',horizontalalignment='center',verticalalignment='center')
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
            disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
            s_disk = s_all[disk]
            cold = f.LowPass('temp', '15000 K') # nur kaltes gas
            s = s_disk.g[cold]
            gas_mass.append(s.g['mass'])
            time.append(s.properties['time'].in_units('Gyr'))

        ax.plot(time,gas_mass,color=color[k],lw=1, linestyle = '-', label = label[k])
        #ax.set_ylim(-0.1,4.45)
        #ax.set_xlim(0.1,1.6)
for n in range(2):
    if n==0: ax.legend(fontsize=12, loc='center',bbox_to_anchor=(-0.5, 0.8))

plt.savefig('gas_mass_time_both.pdf', bbox_inches='tight')
