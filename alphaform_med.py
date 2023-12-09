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
from pynbody import util
import os, struct

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

def starlog(filename):
    f = util.open_(filename, "rb")
    size = struct.unpack(">i", f.read(4))
    iSize = size[0]
    datasize = os.path.getsize(filename) - f.tell()
    datasize % iSize
    file_structure = np.dtype({'names': ("iord", "iorderGas", "tform",
                                                        "x", "y", "z",
                                                        "vx", "vy", "vz",
                                                        "massform", "rhoform", "tempform",
                                                        "alphaform", "epsilonform"),
                                              'formats': ('i4', 'i4', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8', 'f8',
                                                          'f8', 'f8')})
    return np.fromstring(f.read(datasize), dtype=file_structure).byteswap()

def massweight(array_x, array_y, array_key, array_mass, b, range_tuple):
    hist, xbin, ybin = np.histogram2d(array_x,array_y,weights=array_key*array_mass, bins=b, range = range_tuple)
    mass, xbin, ybin = np.histogram2d(array_x,array_y,weights=array_mass, bins=b, range=range_tuple)
    return hist/mass

def load_sim_faceon(mod):
    s = pynbody.load('../med'+'_'+mod+'_iso/' + 'med.01000')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    if (mod == 'federrath_tempcut' or mod == 'federrath_new'):
        filename = '../med_'+mod+'_iso/med.starlog'
        g = starlog(filename)
        key.append(g['alphaform'])
        x.append(g['x'])
        y.append(g['y'])
    else:
        key.append(s.g['alphaform'])
        x.append(s.g['x'])
        y.append(s.g['y'])

        
model = ['semenov', 'padoan', 'evans', 'federrath_new']
for m in model:
    load_sim_faceon(m)


print(len(key))
print(len(x))
print(len(y))

# Titel immer zu bearbeiten
titlelist = [r'a) Semenov et al. (2016)', r'b) Padoan et al. (2012)',  r'c) Evans et al. (2022)', r'd) Federrath & Klessen (2012)' + '\n' + 'without temperature cut', '', '', '', '',]
vmin_list = [0, 0, 0, 0]
vmax_list = [1, 1, 1, 100]
fig = plt.figure(figsize = (10, 7))
gs0 = gd.GridSpec(2, 2)
for n in range(4):
    ax = fig.add_subplot(gs0[n])
    print(len(key[n]))
    print(len(x[n]))
    print(len(y[n]))
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], key[n], statistic='mean', bins=bins, range = ((-50, 50), (-50,50)))
    im = ax.imshow(hist, extent=(-50,50,-50,50), cmap='CMRmap_r', vmin=vmin_list[n], vmax=vmax_list[n], origin='lower')
    ax.set_xlim(-19.99, 19.99)
    ax.set_ylim(-19.99, 19.99)
    ax.text(0.5, 0.88, titlelist[n], horizontalalignment='center', transform=ax.transAxes)
    ax.set_xticklabels([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'virial parameter $\alpha$', size=12)

fig.suptitle('Virial Parameter in SF regions (med resolution)')
plt.savefig('alphaform_med.pdf', bbox_inches='tight')
plt.clf()

