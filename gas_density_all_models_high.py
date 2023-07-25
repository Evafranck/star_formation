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


s = pynbody.load('../high_master_iso/high.01000')
s.physical_units()
s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))

s2 = pynbody.load('../high_padoan_iso/high.01000')
s2.physical_units()
s2.g['n'] = s2.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))

s3 = pynbody.load('../high_semenov_iso/high.01000')
s3.physical_units()
s3.g['n'] = s3.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))

s4 = pynbody.load('../high_master_iso/high.01000')
s4.physical_units()
s4.g['n'] = s4.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))

key = [s.g['n'], s2.g['n'], s3.g['n'], s4.g['n'], s.g['n'], s2.g['n'], s3.g['n'], s4.g['n']]
x = [s.g['x'], s2.g['x'], s3.g['x'], s4.g['x'], s.g['x'], s2.g['x'], s3.g['x'], s4.g['x']]
y = [s.g['y'], s2.g['y'], s3.g['y'], s4.g['y'], s.g['y'], s2.g['y'], s3.g['y'], s4.g['y']]
titlelist = ['a) Threshold-based model', 'b) Padoan et al. (2012)','c) Semenov et al. (2016)', 'd) Evans et al. (2022)', '', '', '', '']
fig = plt.figure(figsize = (12,3.73))
gs0 = gd.GridSpec(2, 4, figure=fig, height_ratios = [1, 0.258], width_ratios = [1, 1, 1, 1.072])
gs0.update(hspace=0.00, wspace=0.00)


for n in range(4):
    ax = fig.add_subplot(gs0[n])
    #print(s2.s['n'].max())
    #print(s2.s['n'].min())
    #print(s2.s['x'].max())
    hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=600, range = ((-50, 50), (-50,50)))
    im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = -1.9, vmax = 5)
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
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-19.99, 19.99)
    ax.set_aspect(1./ax.get_data_ratio())

for n in range(4, 8):
    ax = fig.add_subplot(gs0[n])
    pynbody.analysis.angmom.sideon(s)
    pynbody.analysis.angmom.sideon(s2)
    pynbody.analysis.angmom.sideon(s3)
    pynbody.analysis.angmom.sideon(s4)
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(90)
    hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=600, range = ((-50, 50), (-50,50)))
    im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base, vmin = -2, vmax = 5)
    
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
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-5, 5)

plt.savefig('Gas_iso_com_all_high.pdf', bbox_inches='tight')
plt.clf()

