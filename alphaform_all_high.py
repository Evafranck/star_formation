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
x = []
y = []
mass = []

def massweight(array_x, array_y, array_key, array_mass, b, range_tuple):
    hist, xbin, ybin = np.histogram2d(array_x,array_y,weights=array_key*array_mass, bins=b, range = range_tuple)
    mass, xbin, ybin = np.histogram2d(array_x,array_y,weights=array_mass, bins=b, range=range_tuple)
    return hist/mass

def load_sim_faceon(mod):
    s = pynbody.load('../high'+'_'+mod+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    key.append(s.g['alphaform'])
    print((s.g['alphaform']).max())
    print(np.median(s.g['alphaform']))
    x.append(s.g['x'])
    y.append(s.g['y'])
    mass.append(s.g['mass'])

def load_sim_sideon(mod):
    s = pynbody.load('../high'+'_'+mod+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.sideon(s)
    s.physical_units()
    key.append(s.g['alphaform'])
    x.append(s.g['x'])
    y.append(s.g['y'])
    mass.append(s.g['mass'])


model = ['evans', 'padoan', 'semenov', 'federrath']
for m in model:
    load_sim_faceon(m)
for m in model:    
    load_sim_sideon(m)

# Titel immer zu bearbeiten
titlelist = [r'a) Evans et al. (2022)', r'b) Padoan et al. (2012)', r'c) Semenov et al. (2016)', r'd) Federrath et al. (2014)', '', '', '', '',]

fig = plt.figure(figsize = (12, 3.85))
gs0 = gd.GridSpec(2, 4, height_ratios = [1, 0.3], width_ratios = [1, 1, 1, 1.066])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(8):
    # face-on
    if (n<4):
        ax = fig.add_subplot(gs0[n])
        #hist, xbin, ybin = np.histogram2d(x[n], y[n], weights=key[n], bins=600, range = ((-50, 50), (-50,50)))
        im = ax.imshow(massweight(x[n], y[n], key[n], mass[n], 600, ((-20,20),(-20,20))), extent=(-50,50,-50,50), cmap='CMRmap_r', vmax = 4)
        ax.set_xlim(-19.99, 19.99)
        ax.set_ylim(-19.99, 19.99)
        ax.text(0.5, 0.88, titlelist[n], horizontalalignment='center', transform=ax.transAxes)
        ax.set_xticklabels([])
        if (n == 3):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'virial parameter $\alpha$', size=12)
        if (n == 0):
            ax.set_ylabel('y [kpc]', fontsize = 12)

        else:
            ax.set_yticklabels([])

    # side-on
    if (n>3): 
        ax = fig.add_subplot(gs0[n])
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(90)
        im = ax.imshow(massweight(x[n], y[n], key[n], mass[n], 600, ((-20,20),(-10,10))), extent=(-50,50,-50,50), cmap='CMRmap_r', transform = rot+base, vmax = 4)
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

fig.suptitle('Virial Parameter in SF regions (high resolution)')
plt.savefig('alphaform_all_high.pdf', bbox_inches='tight')
plt.clf()

