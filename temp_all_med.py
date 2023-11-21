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
tempform = []
x = []
y = []
x_s = []
y_s = []

def load_sim_faceon(mod):
    s_all = pynbody.load('../med'+'_'+mod+'_iso/' + 'med.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s = s_all[disk]
    key.append(s.g['temp'])
    #tempform.append(s.s['tempform'])
    print(s.g['temp'].max())
    print(s.g['temp'].min())
    x.append(s.g['x'])
    y.append(s.g['y'])
    x_s.append(s.s['x'])
    y_s.append(s.s['y'])

def load_sim_sideon(mod):
    s_all = pynbody.load('../med'+'_'+mod+'_iso/' + 'med.01000')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s = s_all[disk]
    key.append(s.g['temp'])
    #tempform.append(s.s['tempform'])
    x.append(s.g['x'])
    y.append(s.g['y'])
    x_s.append(s.s['x'])
    y_s.append(s.s['y'])


model = ['master', 'semenov', 'evans', 'federrath']
for m in model:
    load_sim_faceon(m)
for m in model:    
    load_sim_sideon(m)

# Titel immer zu bearbeiten
titlelist = [r'a) Threshold-based model', r'b) Semenov et al. (2016)', r'c) Evans et al. (2022)', r'd) Federrath et al. (2014)', '', '', '', '',]

fig = plt.figure(figsize = (12, 3.85))
gs0 = gd.GridSpec(2, 4, height_ratios = [1, 0.3], width_ratios = [1, 1, 1, 1.07])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(8):
    # face-on
    if (n<4):
        ax = fig.add_subplot(gs0[n])
        hist, xbin, ybin = np.histogram2d(x[n], y[n], weights=key[n], bins=300, range = ((-30, 30), (-30,30)))
        #histform, xbins, ybins = np.histogram2d(x_s[n], y_s[n], weights=tempform[n], bins=400, range = ((-30, 30), (-30,30)))
        im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='seismic', vmin = 2, vmax = 8)
        ax.set_xlim(-19.99, 19.99)
        ax.set_ylim(-19.99, 19.99)
        ax.text(0.5, 0.88, titlelist[n], horizontalalignment='center', transform=ax.transAxes)
        ax.set_xticklabels([])
        if (n == 3):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(T) [K]', size=12)
        if (n == 0):
            ax.set_ylabel('y [kpc]', fontsize = 12)

        else:
            ax.set_yticklabels([])

    # side-on
    if (n>3): 
        ax = fig.add_subplot(gs0[n])
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(90)
        hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=300, range = ((-30, 30), (-30,30)))
        im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='seismic', transform = rot+base, vmin = 2, vmax = 8)
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


fig.suptitle('Gas temperature (med resolution)')
plt.savefig('temp_all_med.pdf', bbox_inches='tight')
plt.clf()

