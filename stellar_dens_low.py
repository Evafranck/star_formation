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


plt.rc('axes', linewidth=1)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 1
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 1
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

young = f.LowPass('age', '1.49 Gyr')
#young = f.LowPass('age', '100 Myr')
bins = 300

density = []
x = []
y = []

# Create a list of simulation paths
#simulations = ['../threshold', '../federrath', 'master', 'federrath']
simulations = ['master', 'padoan', 'federrath', 'federrath_tempcut']


# Load a slice of the simulation snapshots faceon and sideon
def load_sim_faceon(mod):
    if (mod == '../threshold' or mod == '../federrath'):
        s_all = pynbody.load(mod + '/halo.00128')
    else:
        s_all = pynbody.load('../low'+'_'+mod+'_iso/' +'low.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    s_disk = s_all[disk]
    s = s_all.s[young]
    s['n'] = s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s['n'])
    x.append(s['x'])
    y.append(s['y'])
    print(mod, 'dens_min = ', s['n'].min(), 'dens_max = ', s['n'].max(), 'dens_mean = ', np.mean(s['n']))


def load_sim_sideon(mod):
    if (mod == '../threshold' or mod == '../federrath'):
        s_all = pynbody.load(mod + '/halo.00128')
    else:
        s_all = pynbody.load('../low'+'_'+mod+'_iso/' +'low.01000')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    s_disk = s_all[disk]
    s = s_all.s[young]
    s['n'] = s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s['n'])
    x.append(s['x'])
    y.append(s['y'])
    print(mod, 'dens_min = ', s['n'].min(), 'dens_max = ', s['n'].max(), 'dens_mean = ', np.mean(s['n']))

for m in simulations:
    load_sim_faceon(m)
for m in simulations:    
    load_sim_sideon(m)
    
titlelist = ['a) Threshold-based model', 'b) Padoan et al. (2012)', 'c) Federrath & Klessen (2012)' + '\n' + 'with temperature cut', 'd) Federrath & Klessen (2012)' + '\n' + 'without temperature cut', '', '', '', '']
#titlelist = ['Threshold-based model' + '\n' + 'Jakob Herpichs ICs', 'Federrath & Klessen (2012)' + '\n' + 'Jakob Herpichs ICs', 'Threshold-based model' + '\n' + 'AGORA ICs', 'Federrath & Klessen (2012)' + '\n' + 'AGORA ICs', '', '', '', ''] # r'Hopkins et al. (2013)' + '\n' + 'with temperature cut', r'Hopkins et al. (2013)' + '\n' + 'without temperature cut'] # 'Federrath & Klessen (2012)' + '\n' + 'with temperature cut']

fig = plt.figure(figsize = (12,3.73))
gs0 = gd.GridSpec(2, 4, figure=fig, height_ratios = [1, 0.258], width_ratios = [1, 1, 1, 1.072])
gs0.update(hspace=0.00, wspace=0.00)


for n in range(4):
    ax = fig.add_subplot(gs0[n])
    #print(s2.s['n'].max())
    #print(s2.s['n'].min())
    #print(s2.s['x'].max())
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], density[n], statistic='mean', bins=bins, range = ((-50, 50), (-50,50)))
    #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=density[n], bins=600, range = ((-50, 50), (-50,50)))
    im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = 0.01, vmax = 3)
    #im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = 0.1, vmax = 5)
    if n == 3:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(n) [particles $cm^{-3}$]', size=14)
    if n == 0:    
        ax.set_ylabel('y [kpc]', fontsize = 14)
    else:
        ax.set_yticklabels([])
    
    ax.set_title(titlelist[n], fontsize = 11)
    ax.set_xlabel('x [kpc]', fontsize = 14)
    ax.set_xlim(-14.99, 14.99)
    ax.set_ylim(-14.99, 14.99)
    ax.set_aspect(1./ax.get_data_ratio())

for n in range(4, 8):
    ax = fig.add_subplot(gs0[n])
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(90)
    #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=density[n], bins=600, range = ((-50, 50), (-50,50)))
    #im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base, vmin = 0, vmax = 5)
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], density[n], statistic='mean', bins=bins, range = ((-50, 50), (-50,50)))
    im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base, vmin = 0.01, vmax = 3)
    
    if n == 7:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical')
    if n == 4:
        ax.set_ylabel('z [kpc]', fontsize = 14)
    else:
        ax.set_yticklabels([])

    #ax.set_title(titlelist[n], fontsize = 11)
    ax.set_xlabel('x [kpc]', fontsize = 14)
    ax.set_xlim(-14.99, 14.99)
    ax.set_ylim(-3.75, 3.75)

plt.savefig('stellar_dens_low.pdf', bbox_inches='tight')
plt.clf()

