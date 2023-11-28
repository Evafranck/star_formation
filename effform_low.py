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
plt.rcParams['axes.edgecolor'] = 'k'#'gray'
#plt.rcParams['axes.grid'] = True
#plt.rcParams['grid.color'] = 'lightgray'
#plt.rcParams['grid.linestyle'] = 'dashed' #dashes=(5, 1)
#plt.rcParams['lines.dashed_pattern'] = 10, 3
#plt.rcParams['grid.linewidth'] = 0.5
#plt.rcParams['axes.facecolor'] = 'whitesmoke'
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fancybox'] = True
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.shadow'] = False
plt.rcParams['legend.edgecolor'] = 'darkgray'
plt.rcParams['patch.linewidth'] = 0.5

key = []
x = []
y = []
bins = 300

def load_sim_faceon(mod):
    s = pynbody.load('../low'+'_'+mod+'_iso/' + 'low.01000')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    key.append(s.g['effform'])
    print((s.g['effform']).max())
    print(np.median(s.g['effform']))
    x.append(s.g['x'])
    y.append(s.g['y'])

def load_sim_sideon(mod):
    s = pynbody.load('../low'+'_'+mod+'_iso/' + 'low.01000')
    pynbody.analysis.angmom.sideon(s)
    s.physical_units()
    key.append(s.g['effform'])
    x.append(s.g['x'])
    y.append(s.g['y'])


model = ['evans', 'padoan', 'federrath_tempcut', 'federrath_new']
for m in model:
    load_sim_faceon(m)
for m in model:    
    load_sim_sideon(m)

# Titel immer zu bearbeiten
titlelist = [r'a) Threshold-based model', r'b) Padoan et al. (2012)',  r'c) Federrath et al. (2014)' + '\n' + 'with temperature cut', r'd) Federrath et al. (2014)' + '\n' + 'without temperature cut', '', '', '', '',]

fig = plt.figure(figsize = (12, 3.85))
gs0 = gd.GridSpec(2, 4, height_ratios = [1, 0.3], width_ratios = [1, 1, 1, 1.066])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(8):
    # face-on
    if (n<4):
        ax = fig.add_subplot(gs0[n])
        hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], key[n], statistic='mean', bins=bins, range = ((-50, 50), (-50,50)))
        im = ax.imshow(hist, extent=(-50,50,-50,50), cmap='CMRmap_r', vmin= 0,vmax = 1)
        ax.set_xlim(-19.99, 19.99)
        ax.set_ylim(-19.99, 19.99)
        ax.text(0.5, 0.88, titlelist[n], horizontalalignment='center', transform=ax.transAxes)
        ax.set_xticklabels([])
        if (n == 3):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'effiency $\epsilon$', size=12)
        if (n == 0):
            ax.set_ylabel('y [kpc]', fontsize = 12)

        else:
            ax.set_yticklabels([])

    # side-on
    if (n>3): 
        ax = fig.add_subplot(gs0[n])
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(90)
        hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], key[n], statistic='mean', bins=bins, range = ((-50, 50), (-50,50)))
        im = ax.imshow(hist, extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base,  vmin= 0,vmax = 1)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlim(-19.99, 19.99)
        ax.set_ylim(-5.99, 5.99)
        ax.set_xlabel('x [kpc]', fontsize = 12)        
        if (n == 7):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical')

        if (n == 4):    
            ax.set_ylabel('z [kpc]', fontsize = 12)
        
        else:
            ax.set_yticklabels([])

fig.suptitle('SFE in SF regions (low resolution)')
plt.savefig('effform_low.pdf', bbox_inches='tight')
plt.clf()

