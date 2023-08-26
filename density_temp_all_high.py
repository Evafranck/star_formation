import pynbody
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import transforms
import numpy as np
from pynbody.analysis import profile
import matplotlib.gridspec as gd
from pynbody import units as units
from pynbody import array
import pynbody.filt as f


density = []
temp = []
mass = []
massform = []
model = ['master', 'semenov', 'evans', 'federrath']
titlelist = ['Threshold-based model', 'Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']

for n in range(4):
    s = pynbody.load('../high'+'_'+model[n]+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s.g['n'])
    temp.append(s.g['temp'])
    mass.append(s.g['mass'])
    massform.append(s.g['massform'])

fig = plt.figure(figsize = (9.92,10))
gs0 = gd.GridSpec(2, 2, figure=fig, width_ratios=[1,1], height_ratios=[1,1])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(4):
    ax = fig.add_subplot(gs0[n])
    hist, xbin, ybin = np.histogram2d(np.log10(density[n]), np.log10(temp[n]), weights=mass[n], bins=400, range = ((-18, 5), (1,9)))
    histform, xbins, ybins = np.histogram2d(np.log10(density[n]), np.log10(temp[n]), weights=massform[n], bins=400, range = ((-18, 5), (1,9)))
    im = ax.imshow(hist, cmap = 'magma_r', extent=(-18, 5, 1,9), norm = matplotlib.colors.LogNorm())
    #im = ax.imshow(histform, cmap = 'magma_r', extent=(-18, 5, 1,9), norm = matplotlib.colors.LogNorm())
    ax.contour(np.flipud(histform),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                linewidths=0.5, cmap = plt.cm.viridis, levels = [1e5, 1e6])
    ax.set_xlim(-18, 5)
    ax.set_ylim(1.5,7.9)
    ax.vlines(-3, 1.5, 7.9, ls = 'dashed', color = 'grey', linewidths = 0.5) # 0.01% of density threshold (10 particles/cm^3)
    ax.hlines(4, -18, 5, ls = 'dashed', color = 'grey', lw = 0.5) # separates hot and cold gas 
    ax.text(0.5, 0.9, titlelist[n], fontsize = 16, transform = ax.transAxes, horizontalalignment = 'center')
    if (n == 0 or n == 1):
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r'log($n_H$ [cm$^{-3}$])', fontsize = 12)
    if (n == 1 or n == 3):
        ax.set_yticklabels([])
        ax.tick_params(left=False)
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes('right', size = '5%', pad = 0.05)
        #fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'Mass [M$_{\odot}$]', size=12)
    else:
        ax.set_ylabel('log(T [K])', fontsize = 12)
    ax.set_aspect(1./ax.get_data_ratio())
    if (n == 2):
        cax = plt.axes([0.07, 0.1, 0.02, 0.15])
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'Mass [M$_{\odot}$]', size=12)
fig.tight_layout()
plt.savefig('density_temp_all_high.pdf')
