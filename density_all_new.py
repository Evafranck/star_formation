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

key = []
x = []
y = []

def load_sim(mod, res):
    s = pynbody.load('../' + res + '_' + mod + '_iso/' + res + '.01000')
    s.physical_units()
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    key.append(s.g['n'])
    x.append(s.g['x'])
    y.append(s.g['y'])
    pynbody.analysis.angmom.sideon(s)
    key.append(s.g['n'])
    x.append(s.g['x'])
    y.append(s.g['y'])

model = ['master', 'evans', 'padoan', 'semenov']
resolution = ['low', 'med', 'high']
for m in model:
    for r in resolution:
        load_sim(m, r)

# Titel immer zu bearbeiten
titlelist = [r'a) Threshold-based model' + '\n' + 'Low Resolution', r'b) Padoan et al. (2012)'+ '\n' + 'Low Resolution', r'c) Semenov et al. (2016)'+ '\n' + 'Low Resolution', r'd) Evans et al. (2022)'+ '\n' + 'Low Resolution', '', '', '', '',
             r'a) Threshold-based model'+ '\n' + 'Medium Resolution', r'b) Padoan et al. (2012)'+ '\n' + 'Medium Resolution', r'c) Semenov et al. (2016)'+ '\n' + 'Medium Resolution', r'd) Evans et al. (2022)'+ '\n' + 'Medium Resolution', '', '', '', '',
             r'a) Threshold-based model' + '\n' + 'High Resolution', r'b) Padoan et al. (2012)'  + '\n' + 'High Resolution', r'c) Semenov et al. (2016)'  + '\n' + 'High Resolution', r'd) Evans et al. (2022)' + '\n' + 'High Resolution', '', '', '', '',]

fig = plt.figure(figsize = (8, 12))
gs0 = gd.GridSpec(6, 4, figure=fig, height_ratios = [1, 0.258, 1, 0.258, 1, 0.258, 1, 0.258, 1, 0.258], width_ratios = [1, 1, 1, 1, 1, 1.072])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(24):
    # face-on
    if (n<4 or (n>8 and n<12) or (n>15 and n<20)):
        ax = fig.add_subplot(gs0[n])
        hist, xbin, ybin = np.histogram2d(x[n], y[n], weights=key[n], bins=600, range = ((-50, 50), (-50,50)))
        im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = -1.9, vmax = 5)
    
    # side-on
        ax = fig.add_subplot(gs0[n+4])
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(90)
        hist, xbin, ybin = np.histogram2d(x[n+1], y[n+1],weights=key[n+1], bins=600, range = ((-50, 50), (-50,50)))
        im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base, vmin = -2, vmax = 5)
    
    if n == 3:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(n) [particles $cm^{-3}$]', size=14)
        
    if (n == 0 or n == 8 or n == 16):    
        ax.set_ylabel('y [kpc]', fontsize = 14)
        ax.set_yticklabels([])
        
    if (n == 4 or n == 12 or n == 20):    
        ax.set_ylabel('z [kpc]', fontsize = 14)
        ax.set_yticklabels([])
        
    if (n>19):
        ax.set_xlabel('x [kpc]', fontsize = 14)
    
    ax.set_title(titlelist[n], fontsize = 11)
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-19.99, 19.99)
    ax.set_aspect(1./ax.get_data_ratio())

plt.savefig('Gas_iso_com_all_high.pdf', bbox_inches='tight')
plt.clf()

