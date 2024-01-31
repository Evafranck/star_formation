import pynbody
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gd
from matplotlib import transforms
from pynbody import units
import pynbody.filt as f
import glob
import scipy.stats

plt.rc('axes', linewidth=0.5)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 0.5
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.major.size'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.size'] = 0.5
plt.rcParams['ytick.minor.width'] = 0.5
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fancybox'] = True
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.shadow'] = False
plt.rcParams['legend.edgecolor'] = 'darkgray'
plt.rcParams['patch.linewidth'] = 0.5

key = []
x = []
y = []
bins = 150

def load_sim_faceon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s_disk = s_all[disk]
    cold = f.LowPass('temp', '30000 K') # nur kaltes gas
    s = s_disk.g[cold]    
    key.append(np.sqrt(s.g['v2'])) #geschwindigkeit
    x.append(s.g['x'])
    y.append(s.g['y'])

def load_sim_sideon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s_disk = s_all[disk]
    cold = f.LowPass('temp', '30000 K') # nur kaltes gas
    s = s_disk.g[cold]    
    key.append(np.sqrt(s.g['v2'])) # geschwindigkeit
    x.append(s.g['x'])
    y.append(s.g['y'])

model = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan']

titlelist = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r'$\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)']

for m in model:
    load_sim_faceon(m)
for m in model:    
    load_sim_sideon(m)

fig = plt.figure(figsize = (14.9, 3.85))
gs0 = gd.GridSpec(2, 5, height_ratios=[1,0.3], width_ratios=[1,1,1,1,1.07])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(10):
    # face-on
    if (n<5):
        ax = fig.add_subplot(gs0[n])
        hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], key[n], statistic='std', bins=bins, range = ((-30, 30), (-30,30)))
        im = ax.imshow(hist, extent=(-30,30,-30,30), cmap='CMRmap_r', vmin = 0.1, vmax = 55)
        ax.set_xlim(-19.99, 19.99)
        ax.set_ylim(-19.99, 19.99)
        ax.text(0.5, 0.88, titlelist[n], horizontalalignment='center', transform=ax.transAxes)
        ax.set_xticklabels([])
        if (n == 4):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'$\sigma$ [km/s]', size=12)
        if (n == 0):
            ax.set_ylabel('y [kpc]', fontsize = 12)

        else:
            ax.set_yticklabels([])

    # side-on
    if (n>4): 
        ax = fig.add_subplot(gs0[n])
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(90)
        # binned statistics nochmal anschauen
        hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], key[n], statistic='std', bins=bins, range = ((-30, 30), (-30,30)))
        im = ax.imshow(hist, extent=(-30,30,-30,30), cmap='CMRmap_r', transform = rot+base, vmin = 0, vmax = 55)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlim(-19.99, 19.99)
        ax.set_ylim(-5.99, 5.99)
        ax.set_xlabel('x [kpc]', fontsize = 12)        
        if (n == 9):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical')

        if (n == 5):    
            ax.set_ylabel('z [kpc]', fontsize = 12)
        
        else:
            ax.set_yticklabels([])


fig.suptitle('Gas velocity dispersion')
plt.savefig('vdisp_map.pdf', bbox_inches='tight')
plt.clf()