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
from pynbody import filt
import os, struct
from pynbody import util


density = []
temp = []
mass = []
dens_sf = []
temp_sf = []
mass_sf = []
model = ['../low_master_iso', '../low_padoan_iso', '../low_federrath_iso', '../low_federrath_tempcut_iso']
titlelist = ['Threshold-based model', 'Padoan et al. (2012)', 'Federrath & Klessen (2012)' + '\n' + 'without temperature cut', 'Federrath & Klessen (2012)' + '\n' + 'with temperature cut']
#model = ['../threshold', '../federrath', '../low_master_iso', '../low_federrath_new_iso'] # '../low_hopkins_iso', '../low_federrath_tempcut_iso']
#titlelist = ['Threshold-based model' + '\n' + 'Jakob Herpichs ICs', 'Federrath & Klessen (2012)' + '\n' + 'Jakob Herpichs ICs', 'Threshold-based model' + '\n' + 'AGORA ICs', 'Federrath & Klessen (2012)' + '\n' + 'AGORA ICs'] # r'Hopkins et al. (2013)' + '\n' + 'with temperature cut', r'Hopkins et al. (2013)' + '\n' + 'without temperature cut'] # 'Federrath & Klessen (2012)' + '\n' + 'with temperature cut']


def starlog(filename):
    f = util.open_(filename, "rb")
    size = struct.unpack(">i", f.read(4))
    iSize = size[0]
    datasize = os.path.getsize(filename) - f.tell()
    datasize % iSize
    file_structure = np.dtype({'names': ("iord", "iorderGas", "tform",
                                                        "x", "y", "z",
                                                        "vx", "vy", "vz",
                                                        "massform", "rhoform", "tempform",
                                                        "alphaform", "epsilonform"),
                                              'formats': ('i4', 'i4', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8')})
    return np.fromstring(f.read(datasize), dtype=file_structure).byteswap()

def load_sim_faceon(mod):
    if (mod == '../threshold' or mod == '../federrath'):
        s_all = pynbody.load(mod + '/halo.00128')
    else:
        s_all = pynbody.load(mod + '/low.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = filt.LowPass('r', '30 kpc') & filt.BandPass('z', '-5 kpc', '5 kpc')
    s = s_all[disk]
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s.g['n'])
    temp.append(s.g['temp'])
    mass.append(s.g['mass'])
    if (mod == '../low_federrath_tempcut_iso' or mod == '../low_federrath_new_iso' or mod == '../low_hopkins_tempcut_iso' or mod == '../low_hopkins_iso'):
        new = filt.LowPass('age', '1 Gyr')
        filename = mod + '/low.starlog'
        g_tempcut = starlog(filename)
        dens_sf.append(g_tempcut['rhoform']*40.8)
        temp_sf.append(g_tempcut['tempform'])
        mass_sf.append(g_tempcut['massform']*10**9)
    elif (mod == '../federrath'):
        new = filt.LowPass('age', '1 Gyr')
        filename = mod + '/halo.starlog'
        g_tempcut = starlog(filename)
        dens_sf.append(g_tempcut['rhoform']*40.8)
        temp_sf.append(g_tempcut['tempform'])
        mass_sf.append(g_tempcut['massform']*10**9)
    elif (mod == '../low_master_iso' or mod == '../threshold'):
        new = filt.LowPass('age', '1 Gyr')
        s2 = s.s[new]
        s2.s['n_sf'] = s2.s['rhoform'].in_units('kg cm^-3')/(1.673*10**(-27))
        dens_sf.append(s2.s['n_sf'])
        temp_sf.append(s2.s['tempform'])
        mass_sf.append(s2.s['massform'])
    else:
        dens_sf.append([])
        temp_sf.append([])
        mass_sf.append([])
    
for m in model:
    load_sim_faceon(m)
        
print(dens_sf[0])
print(mass_sf[0])  
print(dens_sf[1])
print(mass_sf[1])
print(dens_sf[2])
print(mass_sf[2]) 
print(dens_sf[3])
print(temp_sf[3])
print(mass_sf[3]) 

fig = plt.figure(figsize = (9.92,10))
gs0 = gd.GridSpec(2, 2, figure=fig, width_ratios=[1,1], height_ratios=[1,1])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(4):
    ax = fig.add_subplot(gs0[n])
    hist, xbin, ybin = np.histogram2d(np.log10(density[n]), np.log10(temp[n]), weights=mass[n], bins=400, range = ((-6, 5), (1.5,7.9)))
    im = ax.imshow(np.rot90(hist), cmap = 'magma_r', extent=[xbin[0],xbin[-1],ybin[0],ybin[-1]], norm = matplotlib.colors.LogNorm(vmin = 10**(4.5), vmax = 10**(7.5)))

    if (n<5):
        #print(np.log10(float(dens_sf[-1])))
        histform, xbins, ybins = np.histogram2d(np.log10(dens_sf[n]), np.log10(temp_sf[n]), weights=mass_sf[n], bins=400, range = ((-6, 5), (1.5,7.9)))
        #im = ax.imshow(np.rot90(histform), cmap = 'magma_r', extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], norm = matplotlib.colors.LogNorm())
        level = [1e5, 1e6]
        colors = ['orchid', 'green']
        strs = [r'$10^5 M_{\rm sun}$', r'$10^6 M_{\rm sun}$']
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
plt.savefig('dens_temp_low.pdf')
