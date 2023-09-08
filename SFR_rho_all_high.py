import pynbody
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
from pynbody.analysis import profile
import matplotlib.gridspec as gd
from pynbody import units as units
from pynbody import array
import pynbody.filt as f


model = ['master', 'padoan', 'semenov', 'evans', 'federrath']
labellist = ['Threshold-based model', 'Padoan et al. (2012)', 'Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']
colorlist = ['blue','orange', 'green', 'red', 'purple']

plt.figure(figsize = (10,10))
plt.title('Kennicutt-Schmidt Relation', fontsize = 16)
plt.xlabel(r'$\Sigma_{\rm{gas}}$ [M$_\odot$/pc$^2$]', fontsize = 14)
plt.ylabel(r'$\Sigma_{\rm{SFR}}$ [M$_\odot$/kpc$^2$/yr]', fontsize = 14)
#plt.xlim(0, 1.5)
#plt.ylim(0, 15)
for n in range(5):
    s = pynbody.load('../high' + '_'+ model[n] + '_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    plt.scatter(np.log10(s.g['rho'].in_units('Msol pc**-3')*1.4), np.log10(s.g['sfr'].in_units('Msol kpc**-2 yr**-1')), s = 0.5, color = colorlist[n], label = labellist[n])
    #sfh(s, label = labellist[n], lw = 1, color = colorlist[n])
plt.legend()
plt.savefig('SFH_rho_all_high.pdf', bbox_inches = 'tight')
