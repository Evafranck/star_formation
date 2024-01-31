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

#young = f.LowPass('age', '1.49 Gyr')
#young = f.LowPass('age', '100 Myr')
bins = 300

density = []
x = []
y = []
model = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan', 'threshold']
titlelist = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r'$\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)', 'platzhalter']

def load_sim_faceon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    #s = s_all.s[young]
    s = s_all.s
    s['n'] = s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s['n'])
    x.append(s['x'])
    y.append(s['y'])
    #print(mod, s.g['rho'].min(), s.g['rho'].max(), np.highian(s.g['rho']))
    print(mod, 'dens_min = ', s['n'].min(), 'dens_max = ', s['n'].max(), 'dens_mean = ', np.mean(s['n']))

def load_sim_sideon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    #s = s_all.s[young]
    s = s_all.s
    s['n'] = s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s['n'])
    x.append(s['x'])
    y.append(s['y'])
    #print(mod, s.g['rho'].min(), s.g['rho'].max(), np.median(s.g['rho']))
    print(mod, 'dens_min = ', s['n'].min(), 'dens_max = ', s['n'].max(), 'dens_mean = ', np.mean(s['n']))
    
for m in model:
    load_sim_faceon(m)
for m in model:    
    load_sim_sideon(m)
    
fig = plt.figure(figsize = (14, 2.9))
gs0 = gd.GridSpec(2, 6, figure=fig, height_ratios = [1, 0.258], width_ratios = [1, 1, 1, 1, 1, 1.07])
gs0.update(hspace=0.00, wspace=0.00)


for n in range(6):
    ax = fig.add_subplot(gs0[n])
    #print(s2.s['n'].max())
    #print(s2.s['n'].min())
    #print(s2.s['x'].max())
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], density[n], statistic='mean', bins=bins, range = ((-30, 30), (-30,30)))
    #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=density[n], bins=600, range = ((-50, 50), (-50,50)))
    im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r', vmin = -1.9, vmax = 2)
    #im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = 0.1, vmax = 5)
    if n == 5:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(n) [particles $cm^{-3}$]', size=10)
    
    if n == 0:    
        ax.set_ylabel('y [kpc]', fontsize = 10)
    else:
        ax.set_yticklabels([])
    
    ax.text(0.5, 0.88, titlelist[n], fontsize = 8, horizontalalignment='center', transform=ax.transAxes)
    ax.set_xlabel('x [kpc]', fontsize = 14)
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-19.99, 19.99)
    ax.set_aspect(1./ax.get_data_ratio())

for n in range(6, 12):
    ax = fig.add_subplot(gs0[n])
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(90)
    #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=density[n], bins=600, range = ((-30, 30), (-30,30)))
    #im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r', transform = rot+base, vmin = 0, vmax = 5)
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], density[n], statistic='mean', bins=bins, range = ((-30, 30), (-30,30)))
    im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r', transform = rot+base, vmin = -1.9, vmax = 2)
    
    if n == 11:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical')
    if n == 6:
        ax.set_ylabel('z [kpc]', fontsize = 10)
    else:
        ax.set_yticklabels([])

    #ax.set_title(titlelist[n], fontsize = 11)
    ax.set_xlabel('x [kpc]', fontsize = 10)
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-5, 5)

plt.savefig('stellar_dens.pdf', bbox_inches='tight')
plt.clf()

