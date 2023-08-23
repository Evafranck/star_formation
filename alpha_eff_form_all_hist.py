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
from scipy.optimize import curve_fit

alpha_evans = []
alpha_padoan = []
alpha_semenov = []
alpha_federrath = []

eff_evans = []
eff_padoan = []
eff_semenov = []
eff_federrath = []

eff_list = [eff_evans, eff_padoan, eff_semenov, eff_federrath]
alpha_list = [alpha_evans, alpha_padoan, alpha_semenov, alpha_federrath]
keylist = [alpha_list, eff_list]
x_list = [(-0.1, 1.1), (-0.1, 1.1)]
label_list = ['Evans et al. (2022)', 'Padoan et al. (2012)', 'Semenov et al. (2016)', 'Federrath et al. (2014)']
title_list = ['Virial Parameter in SF regions (high resolution)', 'SFE in SF regions (high resolution)']
colorlist = ['blue', 'red', 'green', 'purple']
xlabel_list = [r'$\alpha_{\rm{sf}}$', r'$\epsilon_{\rm{sf}}$']
model = ['evans', 'padoan', 'semenov', 'federrath']

for i in range(4):
    s = pynbody.load('../high'+'_'+model[i]+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    alpha_list[i].append(s.g['alphaform'])
    eff_list[i].append(s.g['effform'])
    print((s.g['effform']).max())
    print(np.median(s.g['effform']))
    
fig = plt.figure(figsize = (10,5))
gs0 = gd.GridSpec(1, 2, figure=fig, wspace = 0.4)

for n in range(2):
    ax = fig.add_subplot(gs0[n])
    for i in range(4):
        hist, bins, edges = ax.hist(keylist[n][i], bins = 100, color = colorlist[i], range = x_list[n], label = label_list[i], histtype = 'step', density = True)
    ax.set_aspect(1./ax.get_data_ratio())
    ax.grid(linewidth = 0.3, color = 'grey')
    ax.set_xlim(x_list[n])
    ax.set_title(title_list[n])
    ax.set_xlabel(xlabel_list[n])
plt.legend()

'''        if (n==0):
            mean = np.mean(eff_list[i])
            sigma = np.sqrt(np.var(eff_list[i]))
            print(mean)
            print(sigma)
            ax.plot(bins[np.argmax(hist)+1], hist.max(), "x", color = colorlist[i])
            ax.vlines(x=bins[np.argmax(hist)+1], ymax = hist.max(), ymin = 0, color = colorlist[i], linestyle = 'dashed', linewidth = 1)
    ax.legend(fontsize = 7, loc = loc[n])
    ax.set_title(title[n], wrap = True, fontsize = 12)
    ax.set_xlabel(unit[n])
    ax.set_xlim(range_list[n])
    ax.set_ylim(y_list[n])'''

fig.savefig('alpha_eff_form_high.pdf')
