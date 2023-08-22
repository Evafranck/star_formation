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

#s = pynbody.load('/home/hd/hd_hd/hd_fd249/simulation/agora/low_padoan_cos/g1.37e11.2x2.std')
#s.physical_units()
s = pynbody.load('../low_padoan_iso/LOW.bin')
s2 = pynbody.load('../med_padoan_iso/MED.bin')
s3 = pynbody.load('../high_padoan_iso/HI.bin')
s.physical_units()
s2.physical_units()
s3.physical_units()

s.s['n'] = s.s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
s2.s['n'] = s2.s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
s3.s['n'] = s3.s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))

#key = [s.g['n'], s2.g['n']]
#x = [s.g['x'], s2.g['x']]
#y = [s.g['y'], s2.g['y']]
titlelist = ['a) Low Resolved Isolated Simulation','b) Medium Resolution Isolated Simulation', 'High Resolution Isolated Simulation','','','','Initial star density','Initial gas density', 'Initial star density']
fig = plt.figure(figsize = (18,10))
gs0 = gd.GridSpec(2, 3, figure=fig, height_ratios = [1, 0.25], width_ratios = [1, 1 ,1.058])
gs0.update(hspace=0.00, wspace=0.00)
    
for n in range(3): 
    ax = fig.add_subplot(gs0[n])
    #print(s2.s['n'].max())
    #print(s2.s['n'].min())
    #print(s2.s['x'].max())
    hist1, xbin, ybin = np.histogram2d(s.s['x'], s.s['y'],weights=s.s['n'], bins=400, range = ((-50, 50), (-50,50)))
    hist2, xbin, ybin = np.histogram2d(s2.s['x'], s2.s['y'],weights=s2.s['n'], bins=400, range = ((-50, 50), (-50,50)))
    hist3, xbin, ybin = np.histogram2d(s3.s['x'], s3.s['y'],weights=s3.s['n'], bins=400, range = ((-50, 50), (-50,50)))
    hist = [hist1, hist2, hist3]
    im = ax.imshow(np.log10(hist[n]), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = -1.9, vmax = 6)
    if n == 0:
        ax.set_ylabel('y [kpc]', fontsize = 14)
    else:
        ax.set_yticklabels([])
    if n == 1:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(n) [particles $cm^{-3}$]', size=14)
#hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=100, range = ((-50, 50), (-50,50)
 #   hist3, xbin, ybin = np.histogram2d(s3.s['x'], s3.s['y'],weights=s3.s['n'], bins=400, range = ((-50, 50), (-50,50)))
    ax.set_title(titlelist[n], fontsize = 13)
    ax.set_xlabel('x [kpc]', fontsize = 14)
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-19.99, 19.99)
    ax.set_aspect(1./ax.get_data_ratio())

for n in range (3, 6):
    ax = fig.add_subplot(gs0[n])
    pynbody.analysis.angmom.sideon(s)
    pynbody.analysis.angmom.sideon(s2)
    pynbody.analysis.angmom.sideon(s3)
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(90)
    hist1, xbin, ybin = np.histogram2d(s.s['x'], s.s['y'],weights=s.s['n'], bins=400, range = ((-50, 50), (-50,50)))
    hist2, xbin, ybin = np.histogram2d(s2.s['x'], s2.s['y'],weights=s2.s['n'], bins=400, range = ((-50, 50), (-50,50)))
    hist3, xbin, ybin = np.histogram2d(s3.s['x'], s3.s['y'],weights=s3.s['n'], bins=400, range = ((-50, 50), (-50,50)))
    hist = [hist1, hist2, hist3]
    im = ax.imshow(np.log10(hist[n-3]), extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base, vmin = -2, vmax = 6)
    if n == 3:
        ax.set_ylabel('y [kpc]', fontsize = 14)
    else:
        ax.set_yticklabels([])
    if n == 4:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical')    
    ax.set_title(titlelist[n], fontsize = 11)
    ax.set_xlabel('x [kpc]', fontsize = 14)
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-5, 5)
plt.savefig('density_IC.pdf', bbox_inches='tight')
plt.clf()

