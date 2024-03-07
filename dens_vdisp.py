import pynbody
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter
from matplotlib import transforms
import numpy as np
#from pynbody.analysis import profile
import matplotlib.gridspec as gd
from pynbody import units as units
from pynbody import array
from pynbody import filt
import os, struct
from pynbody import util
import scipy.stats

density = []
eff = []
vdisp = []
temp = []
temperature = []
dens = []
bins = 150
#titlelist = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r' $\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)', 'Platzhalter']

#model = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
#model = ['threshold_alpha008', ] 
model = ['federrath_1e6_alpha008', 'federrath_alpha008', 'federrath_cstar_cut', 'semenov_alpha008','semenov_1e6_alpha008', 'semenov_cstar_cut']
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
    print('units', s_all.s['rho'].units)
    print('units', s_all.g['v_disp'].units)
    pynbody.analysis.angmom.faceon(s_all)
    print(s_all.properties['a'])
    print(s_all.properties['z'])
    print(s_all.properties['time'])
    new = filt.LowPass('age', '15 Myr')
    s_all.physical_units()
    disk = filt.LowPass('r', '20 kpc') & filt.BandPass('z', '-1 kpc', '1 kpc')
    cold = filt.LowPass('temp', '30000 K')
    #s = s_all
    s_disk = s_all[disk]
    s = s_disk.g[cold]
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s.g['n'])
    vdisp.append(s.g['v_disp'])
    temperature.append(s.g['temp'])
    new = filt.LowPass('age', '15 Myr')
    filename = '../' + mod + '/halo.starlog'
    g_tempcut = starlog(filename)
    dens.append(g_tempcut['rhoform']/106.2939218) #*40.8)
    eff.append(g_tempcut['epsilonform'])
    temp.append(g_tempcut['tempform'])
    #print('Efficiency min und max', np.min(g_tempcut['epsilonform']), np.max(g_tempcut['epsilonform']), np.mean(g_tempcut['epsilonform']))
    #print('Temperature min und max', np.min(g_tempcut['tempform']), np.max(g_tempcut['tempform']), np.mean(g_tempcut['tempform']))
    #print('Density min und max', np.min(g_tempcut['rhoform']/106.2939218), np.max(g_tempcut['rhoform']/106.2939218), np.mean(g_tempcut['rhoform']/106.2939218)) 
    print('Velocity Dispersion min und max', np.min(s.g['v_disp']), np.max(s.g['v_disp']), np.mean(s.g['v_disp']))
    print('Temperature min und max', np.min(s.g['temp']), np.max(s.g['temp']), np.mean(s.g['temp']))
    print('Density min und max', np.min(s.g['n']), np.max(s.g['n']), np.mean(s.g['n'])) 
    
for m in model:
    load_sim_faceon(m)

fig = plt.figure(figsize = (18,3))
gs0 = gd.GridSpec(1, 6, figure=fig)
#gs0.update(hspace=0.00, wspace=0.0)

#fig2 = plt.figure()
#axx = plt.subplot(111)

for n in range(6):
    ax = fig.add_subplot(gs0[n])
    #hist, xbin, ybin = np.histogram2d(np.log10(density[n]), vdisp[n], weights=temperature[n], bins=bins)#, range = ((-3, 10), (0, 0.2)))
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(np.log10(density[n]), np.log10(vdisp[n]), temperature[n], statistic='mean', bins=bins, range = ((-3, 3), (0.5, 3)))
    im = ax.imshow(hist.T, origin='lower',cmap = 'coolwarm', extent=[xbin[0],xbin[-1],ybin[0],ybin[-1]], norm = matplotlib.colors.LogNorm(vmin = 10**(2), vmax = 10**6))
    ax.set_xlim(-3, 2.5)
    ax.set_ylim(0.5,3)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y))) 
    ax.text(0.5, 0.9, titlelist[n], fontsize = 10, transform = ax.transAxes, horizontalalignment = 'center')
    ax.set_xlabel(r'log($n_H$ [cm$^{-3}$])', fontsize = 10)
    #if (n != 0):
    #    ax.set_yticklabels([])
    #    ax.tick_params(left=False)
    #else:
    ax.set_ylabel('Velocity Dispersion', fontsize = 10)
    ax.set_aspect(1./ax.get_data_ratio())
    if n == 2:
        cax = plt.axes([0.12, 0.22, 0.01, 0.3])
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'Temperature [K]', size=8)
        #ax.legend(loc = 'lower left', fontsize = 6)
fig.tight_layout()
plt.savefig('dens_vdisp.pdf')
#fig2.legend()
#fig2.savefig('dens_hist.pdf', bbox_inches='tight')