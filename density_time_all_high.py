import pynbody 
import matplotlib.pyplot as plt
import numpy as np
from pynbody.analysis import profile
import matplotlib.gridspec as gd
from pynbody import units as units
from pynbody import array
import pynbody.filt as f

density = []
temp = []
model = ['master', 'padoan', 'semenov', 'evans', 'federrath']
titlelist = ['Threshold-based model', 'Padoan et al. (2012)', 'Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']

for n in model:
    s = pynbody.load('../high'+'_'+model+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s.g['n'])
    temp.append(s.g['temp'])

fig = plt.figure(figsize = (17,5))
gs0 = gd.GridSpec(1, 1, figure=fig)

for n in range(1):
    ax = fig.add_subplot(gs0[n])
    ax.hist2d(np.log10(density[n]), np.log10(temp[n]), bins = 100, cmap = 'viridis', norm = matplotlib.colors.LogNorm())
    ax.set_aspect(1./ax.get_data_ratio())
    ax.set_text(0.5, 0.8, titlelist[n], fontsize = 14, transform = ax.transAxes, horizontalalignment = 'center')
    ax.set_xlabel(r'log($n_H$) [cm$^{-3}$]')
    ax.set_ylabel('log(T) [K]')

plt.colorbar()
plt.savefig('density_temp_all_high.pdf')