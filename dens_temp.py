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
bins = 150
#titlelist = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r' $\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)', 'Platzhalter']

#model = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
#model = ['threshold_alpha008', 'semenov_alpha008','semenov_1e6_alpha008', 'semenov_cstar_cut'] 
model = ['threshold_1e6_alpha008', 'federrath_1e6_alpha008', 'federrath_alpha008', 'federrath_cstar_cut']
titlelist = model

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

def starlog2(filename):
    f = util.open_(filename, "rb")
    size = struct.unpack(">i", f.read(4))
    iSize = size[0]
    datasize = os.path.getsize(filename) - f.tell()
    datasize % iSize
    file_structure = np.dtype({'names': ("iord", "iorderGas", "tform",
                                                        "x", "y", "z",
                                                        "vx", "vy", "vz",
                                                        "massform", "rhoform", "tempform"),
                                              'formats': ('i4', 'i4', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8', 'f8')})
    return np.fromstring(f.read(datasize), dtype=file_structure).byteswap()

def load_sim_faceon(mod):
    s_all = pynbody.load('../' + mod+'/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    print(s_all.properties['a'])
    print(s_all.properties['z'])
    print(s_all.properties['time'])
    new = filt.LowPass('age', '15 Myr')
    if (mod == 'threshold' or mod == 'threshold_alpha008' or mod == 'threshold_1e6_alpha008'):
        print(s_all.s[new]['rhoform'])
        print(s_all.s['rho'].units)
    s_all.physical_units()
    #disk = filt.LowPass('r', '30 kpc') & filt.BandPass('z', '-5 kpc', '5 kpc')
    s = s_all#[disk]
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s.g['n'])
    temp.append(s.g['temp'])
    mass.append(s.g['mass']) 
    if (mod == 'threshold' or mod == 'threshold_alpha008' or mod == 'threshold_1e6_alpha008'):
        new = filt.LowPass('age', '15 Myr')
        filename =  '../' + mod + '/halo.starlog'
        g_tempcut = starlog2(filename)
        s2 = s.s[new]
        s2.s['n_sf'] = s2.s['rhoform'].in_units('kg cm^-3')/(1.673*10**(-27))
        dens_sf.append(s2.s['n_sf'])
        temp_sf.append(s2.s['tempform'])
        mass_sf.append(s2.s['massform'])
        print(mod, 'dens_min = ', s2.s['n_sf'].min(), 'dens_max = ', s2.s['n_sf'].max(), 'dens_mean = ', np.mean(s2.s['n_sf']))
        print('maximum mass = ', s2.s['massform'].max(), 'minimum mass = ', s2.s['massform'].min())

    else:
        new = filt.LowPass('age', '15 Myr')
        filename = '../' + mod + '/halo.starlog'
        g_tempcut = starlog(filename)
        dens_sf.append(g_tempcut['rhoform']/106.2939218) #*40.8)
        temp_sf.append(g_tempcut['tempform'])
        mass_sf.append(g_tempcut['massform']*10**9)
        print(mod, 'dens_min = ', g_tempcut['rhoform'].min()/106.2939218, 'dens_max = ', g_tempcut['rhoform'].max()/106.2939218, 'dens_mean = ', np.mean(g_tempcut['rhoform']/106.2939218))
        print('maximum mass = ', g_tempcut['massform'].max()*10**9, 'minimum mass = ', g_tempcut['massform'].min()*10**9)


for m in model:
    load_sim_faceon(m)
    

fig = plt.figure(figsize = (10,3))
gs0 = gd.GridSpec(1, 4, figure=fig)
gs0.update(hspace=0.00, wspace=0.00)

#fig2 = plt.figure()
#axx = plt.subplot(111)

for n in range(4):
    ax = fig.add_subplot(gs0[n])
    hist, xbin, ybin = np.histogram2d(np.log10(density[n]), np.log10(temp[n]), weights=mass[n], bins=bins, range = ((-6, 8), (1.5,7.9)))
    im = ax.imshow(hist.T, origin='lower',cmap = 'magma_r', extent=[xbin[0],xbin[-1],ybin[0],ybin[-1]], norm = matplotlib.colors.LogNorm(vmin = 10**(3.5), vmax = 10**(8)))

    histform, xbins, ybins = np.histogram2d(np.log10(dens_sf[n]), np.log10(temp_sf[n]), weights=mass_sf[n], bins=bins, range = ((-6, 8), (1.5,7.9)))
    #im = ax.imshow(np.rot90(histform), cmap = 'magma_r', extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], norm = matplotlib.colors.LogNorm())
    level = [1e4, 1e10, 1e11, 5e11, 1e12]
    colors = ['green', 'yellow', 'pink', 'blue', 'red']
    strs = [ r'$10^4 M_{\rm sun}$',  r'$10^{10} M_{\rm sun}$',  r'$10^{11} M_{\rm sun}$', r'$5 \cdot 10^{11} M_{\rm sun}$', r'$10^{12} M_{\rm sun}$']
    cont = ax.contour(histform.T,origin='lower',extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], colors=('r', 'green', 'blue', 'orange', '#afeeee'), linewidths=0.5, levels = level)
    handles, labels = cont.legend_elements()
    ax.set_xlim(-6, 8)
    ax.set_ylim(1.5,7.9)
    ax.vlines(-1.8, 1.5, 7.9, ls = '-', color = 'grey', linewidths = 0.5, label = r'0.01% of $n_{\rm th}$') # 0.01% of density threshold (10 particles/cm^3)
    ax.vlines(1, 1.5, 7.9, ls = '--', color = 'grey', linewidths = 1, label = r'$n_{\rm th}$') # density threshold (10 particles/cm^3)
    ax.hlines(3.5, -6, 8, ls = ':', color = 'grey', lw = 0.5, label = r'T = $3.5 \cdot 10^4$ K') # separates hot and cold gas 
    ax.text(0.5, 0.9, titlelist[n], fontsize = 10, transform = ax.transAxes, horizontalalignment = 'center')
    ax.set_xlabel(r'log($n_H$ [cm$^{-3}$])', fontsize = 10)
    if (n != 0):
        ax.set_yticklabels([])
        ax.tick_params(left=False)
    else:
        ax.set_ylabel('log(T [K])', fontsize = 10)
    ax.set_aspect(1./ax.get_data_ratio())
    if n == 0:
        ax.legend(handles, [r'$10^4 M_{\rm sun}$',  r'$10^{10} M_{\rm sun}$',  r'$10^{11} M_{\rm sun}$', r'$5 \cdot 10^{11} M_{\rm sun}$', r'$10^{12} M_{\rm sun}$'], loc = 'lower right', fontsize = 6)
    if n == 2:
        cax = plt.axes([0.06, 0.2, 0.01, 0.2])
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'Mass [M$_{\odot}$]', size=8)
        ax.legend(loc = 'lower left', fontsize = 6)
fig.tight_layout()
plt.savefig('dens_temp.pdf')
#fig2.legend()
#fig2.savefig('dens_hist.pdf', bbox_inches='tight')