import pynbody
import pylab
import matplotlib.pyplot as plt
import numpy as np
from pynbody import units
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
from scipy import stats
from scipy import constants
import pynbody.filt as f
import matplotlib.gridspec as gd
import pickle
from matplotlib import colorbar
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy.signal import argrelextrema
from pynbody import util
import os, struct

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

efficiency = []
alpha = []
bins = 300


def load_sim_faceon(mod):
    filename = '../'+mod+'/halo.starlog'
    g = starlog(filename)
    s = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units() 
    efficiency.append(g['epsilonform'])
    alpha.append(g['alphaform'])
    #efficiency.append(s['effform'])
    #alpha.append(s['alphaform'])
    print(len(alpha))
    print(len(efficiency))

model = ['federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan']
for m in model:
    load_sim_faceon(m)
    
titlelist = ['Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with $\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)']

range_list = [(4.26, 5.45), (-2, 1), (1, 6)]
y_list = [(0, 8.5), (0, 3.4),(0, 1.1)]
fig = plt.figure(figsize = (20,10))
gs0 = gd.GridSpec(2, 4, figure=fig, hspace=0.2)


for n in range(4):
    ax = fig.add_subplot(gs0[n])
    if n < 6:
        hist, bins, edges = ax.hist(alpha[n], bins = 100, histtype = 'step', density = True)
        ax.set_xlabel(r'$\alpha$', fontsize = 15)
        #ax.set_xlim(0, 300)
    else:
        hist, bins, edges = ax.hist(efficiency[n-4], bins = 100, histtype = 'step', density = True)
        ax.set_xlabel(r'$\epsilon_{\mathrm{ff}}$', fontsize = 15)
        #ax.set_xlim(0, 0.21)
    ax.set_title(titlelist[n], wrap = True, fontsize = 15)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    #ax.set_ylim(y_list[n])
    ax.set_aspect(1./ax.get_data_ratio())
    ax.grid(ls = '--', lw = 0.1, color = 'grey')
fig.tight_layout() 
fig.savefig('alpha_efficiency_histogram.pdf', bbox_inches='tight')