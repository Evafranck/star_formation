import pynbody
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import transforms
import numpy as np
#from pynbody.analysis import profile
import matplotlib.gridspec as gd
from pynbody import units as units
from pynbody import array
import pynbody.filt as f

density = []
temp = []
mass = []
dens_sf = []
temp_sf = []
mass_sf = []
model = ['master', 'semenov', 'evans', 'federrath']
titlelist = ['Threshold-based model', 'Semenov et al. (2016)', 'Evans et al. (2022)', 'Federrath et al. (2014)']

for n in range(4):
    s_all = pynbody.load('../high'+'_'+model[n]+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s = s_all[disk]
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s.g['n'])
    temp.append(s.g['temp'])
    mass.append(s.g['mass'])
    
    if (n==0):
        new = f.LowPass('age', '1 Gyr')
        s2 = s.s[new]
        s2.s['n_sf'] = s2.s['rhoform'].in_units('kg cm^-3')/(1.673*10**(-27))
        dens_sf.append(s2.s['n_sf'])
        temp_sf.append(s2.s['tempform'])
        mass_sf.append(s2.s['massform'])
        
        print('mass_max = ',np.max(mass))
        print('mass_min = ', np.min(mass))
        print('temp_max = ', np.max(temp))
        print('temp_min = ', np.min(temp))
        print('dens_max = ', np.max(density))
        print('dens_min = ', np.min(density))
        print('mass_sf_max = ', np.max(mass_sf))
        print('mass_sf_min = ', np.min(mass_sf))
        print('temp_sf_max = ', np.max(temp_sf))
        print('temp_sf_min = ', np.min(temp_sf))
        print('dens_sf_max = ', np.max(dens_sf))
        print('dens_sf_min = ', np.min(dens_sf))

fig = plt.figure(figsize = (9.92,10))
gs0 = gd.GridSpec(2, 2, figure=fig, width_ratios=[1,1], height_ratios=[1,1])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(4):
    ax = fig.add_subplot(gs0[n])
    hist, xbin, ybin = np.histogram2d(np.log10(density[n]), np.log10(temp[n]), weights=mass[n], bins=300, range = ((-6, 5), (1.5,7.9)))
    im = ax.imshow(np.rot90(hist), cmap = 'magma_r', extent=[xbin[0],xbin[-1],ybin[0],ybin[-1]], norm = matplotlib.colors.LogNorm(vmin = 10**(3.5), vmax = 10**(8)))

    if (n==0):
        histform, xbins, ybins = np.histogram2d(np.log10(dens_sf[n]), np.log10(temp_sf[n]), weights=mass_sf[n], bins=400, range = ((-6, 5), (1.5,7.9)))
        #im = ax.imshow(np.rot90(histform), cmap = 'magma_r', extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], norm = matplotlib.colors.LogNorm())
        level = [1e4, 1e6]
        colors = ['orchid', 'green']
        strs = [r'$10^4 M_{\rm sun}$', r'$10^6 M_{\rm sun}$']
        cont = ax.contour(np.flipud(np.rot90(histform)),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], linewidths=0.5, cmap = plt.cm.PiYG, levels = level)
        list = []
        for i, level in enumerate(level):
            list.append(plt.Line2D([0], [0], ls = '-', lw = 1, color = colors[i]))
        ax.legend(list, strs, loc = 'lower left')
    ax.set_xlim(-6, 5)
    ax.set_ylim(1.5,7.9)
    ax.vlines(-1.8, 1.5, 7.9, ls = '-', color = 'grey', linewidths = 0.5, label = r'0.01% of $n_{\rm th}$') # 0.01% of density threshold (10 particles/cm^3)
    ax.vlines(1, 1.5, 7.9, ls = '--', color = 'grey', linewidths = 1, label = r'$n_{\rm th}$') # density threshold (10 particles/cm^3)
    ax.hlines(3.5, -6, 5, ls = ':', color = 'grey', lw = 0.5, label = r'T = $3.5 \cdot 10^4$ K') # separates hot and cold gas 
    ax.text(0.5, 0.9, titlelist[n], fontsize = 16, transform = ax.transAxes, horizontalalignment = 'center')
    if (n == 0 or n == 1):
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r'log($n_H$ [cm$^{-3}$])', fontsize = 12)
    if (n == 1 or n == 3):
        ax.set_yticklabels([])
        ax.tick_params(left=False)
    else:
        ax.set_ylabel('log(T [K])', fontsize = 12)
    ax.set_aspect(1./ax.get_data_ratio())
    if (n == 2):
        cax = plt.axes([0.07, 0.1, 0.02, 0.15])
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'Mass [M$_{\odot}$]', size=12)
        ax.legend(loc = 'lower right')
fig.tight_layout()
plt.savefig('dens_temp_high.pdf')
