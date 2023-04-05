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

young100 = f.LowPass('age', '100 Myr')

s_ = pynbody.load('../med_master_iso/med.01000')
s2_ = pynbody.load('../med_padoan_iso/med.01000')
s_.physical_units()
s2_.physical_units()
s = s_.s[young100]
s2 = s2_.s[young100]
s.s['n'] = s.s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
s2.s['n'] = s2.s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))

s3_ = pynbody.load('../med_semenov_iso/med.01000')
s4_ = pynbody.load('../med_evans_iso/med.01000')
s3_.physical_units()
s4_.physical_units()
s3 = s3_.s[young100]
s4 = s4_.s[young100]
s3.s['n'] = s3.s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
s4.s['n'] = s4.s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))

key = [s.s['n'], s2.s['n'], s3.s['n'], s4.s['n'], s.s['n'], s2.s['n'], s3.s['n'], s4.s['n']]
x = [s.s['x'], s2.s['x'], s3.s['x'], s4.s['x'], s.s['x'], s2.s['x'], s3.s['x'], s4.s['x']]
y = [s.s['y'], s2.s['y'], s3.s['y'], s4.s['y'], s.s['y'], s2.s['y'], s3.s['y'], s4.s['y']]
titlelist = ['a) Threshold-based model'+'\n'+'stellar age < 100 Myr', 'b) Padoan et al. (2012)' + '\n' + 'stellar age < 100 Myr','c) Semenov et al. (2016)'+'\n'+ 'stellar age < 100 Myr', 'd) Evans et al. (2022)'+'\n'+'stellar age < 100 Myr', '', '', '', '']
fig = plt.figure(figsize = (12,3.73))
gs0 = gd.GridSpec(2, 4, figure=fig, height_ratios = [1, 0.258], width_ratios = [1, 1, 1, 1.072])
gs0.update(hspace=0.00, wspace=0.00)


for n in range(4):
    ax = fig.add_subplot(gs0[n])
    #print(s2.s['n'].max())
    #print(s2.s['n'].min())
    #print(s2.s['x'].max())
    hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=600, range = ((-50, 50), (-50,50)))
    im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = 0.1, vmax = 5)
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
    im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base, vmin = 0, vmax = 5)
    
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

plt.savefig('young100_stars_iso_com_all_med.pdf', bbox_inches='tight')
plt.clf()

