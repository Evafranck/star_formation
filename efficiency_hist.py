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
from matplotlib.ticker import ScalarFormatter
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

eff = []
bins = 200


def load_sim_faceon(mod):
    filename = '../'+mod+'/halo.starlog'
    g = starlog(filename)
    s = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units() 
    eff.append(g['epsilonform'])
    

#model = ['federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan']
#model = ['federrath_1e6_alpha008', 'federrath_alpha008', 'federrath_cstar_cut', 'hopkins_alpha_alpha008', 'hopkins_alpha_padoan', 'hopkins_alpha_padoan_alpha008', 'hopkins_alpha008', 'semenov_1e6_alpha008', 'semenov_alpha008', 'semenov_cstar_cut']
model = ['federrath_1e6_alpha008', 'federrath_alpha008', 'federrath_cstar_cut', 'semenov_1e6_alpha008', 'semenov_alpha008', 'semenov_cstar_cut']

titlelist = model
for m in model:
    load_sim_faceon(m)
#titlelist = ['Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with $\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)']

range_list = [(4.26, 5.45), (-2, 1), (1, 6)]
y_list = [(0, 8.5), (0, 3.4),(0, 1.1)]

fig = plt.figure(figsize = (17,3))
gs0 = gd.GridSpec(1, 6, figure=fig)
#gs0.update(hspace=0.00, wspace=0.00)

for n in range(6):
    ax = fig.add_subplot(gs0[n])
    hist, bins, edges = ax.hist(eff[n], bins = bins, histtype = 'step', density = True)
    ax.set_title(titlelist[n], wrap = True)
    #if n>0 and n!=5:
        #ax.set_yticklabels([])
    #ax.set_xscale('log')
    ax.set_xlabel(r'$\epsilon_{\mathrm{ff}}$')
    #ax.set_xlim(0.01, 5000)
    #ax.set_yscale('log')
    #ax.set_ylim(0, 0.001)
    #plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    #plt.gca().tick_params(axis='y', which='both', bottom=True)
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #ax.tick_params(axis='x', la14belsize=)
    #ax.tick_params(axis='y', labelsize=14)
    #ax.set_ylim(0, 50)
    ax.set_aspect(1./ax.get_data_ratio())
    ax.grid(ls = '--', lw = 0.1, color = 'grey')
fig.tight_layout() 
fig.savefig('eff_histogram.pdf', bbox_inches='tight')