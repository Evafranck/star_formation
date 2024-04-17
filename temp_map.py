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
model_name = 'federrath_padoan' # options: 'hopkins', 'semenov', 'federrath_padoan', '1e6'
x_y_lim = 14.99
z_lim = x_y_lim * 0.266

def load_sim_faceon(mod):
    s_all = pynbody.load('../' +mod + '/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s_disk = s_all[disk]
    cold = f.LowPass('temp', '30000 K')
    s = s_disk.g[cold]
    key.append(s.g['temp'])
    #tempform.append(s.s['tempform'])
    print(s.g['temp'].max())
    print(s.g['temp'].min())
    x.append(s.g['x'])
    y.append(s.g['y'])
    x_s.append(s.s['x'])
    y_s.append(s.s['y'])

def load_sim_sideon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s_disk = s_all[disk]
    cold = f.LowPass('temp', '30000 K')
    s = s_disk.g[cold]
    key.append(s.g['temp'])
    #tempform.append(s.s['tempform'])
    x.append(s.g['x'])
    y.append(s.g['y'])
    x_s.append(s.s['x'])
    y_s.append(s.s['y'])

if model_name == '1e6':
    #models with 10**6 resolution
    model = ['threshold_1e6_alpha008', 'padoan_1e6_alpha008', 'federrath_1e6_alpha008', 'federrath_cstar_cut_1e6', 'semenov_1e6_alpha008', 'semenov_cstar_cut_1e6', 'evans_1e6_alpha008']
    bins = 250
elif model_name == 'hopkins':
    # models that are based on Hopkins et al. (2013)
    model = ['threshold', 'hopkins', 'hopkins_alpha', 'hopkins_alpha008', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
    bins = 200
elif model_name == 'semenov':
    # models that are based on Semenov et al. (2016)
    model = ['threshold_alpha008', 'semenov_alpha008', 'semenov_cstar_cut', 'semenov_alpha008_tcr', 'semenov_alpha008_converging_flow', 'semenov_alpha008_tcool_cut', 'semenov_alpha008_tcool_cut_converging_flow']
    bins = 230
elif model_name == 'federrath_padoan':
    # models that are based on Federrath & Klessen (2012) and Padoan et al. (2012)
    model = ['threshold_alpha008', 'federrath_alpha008', 'federrath_cstar_cut', 'padoan_alpha008']
    bins = 260
    
#titlelist = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r'$\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)', 'platzhalter']
titlelist = model

for m in model:
    load_sim_faceon(m)
for m in model:    
    load_sim_sideon(m)

fig = plt.figure(figsize = (len(model)*2.33, 2.95))
width = np.ones(len(model))
width[-1] = 1.077
gs0 = gd.GridSpec(2, len(model), figure=fig, height_ratios = [1, 0.27], width_ratios = width)
gs0.update(hspace=0.00, wspace=0.00)

for n in range(2*len(model)):
    # face-on
    if (n<len(model)):
        ax = fig.add_subplot(gs0[n])
        hist, xbin, ybin = np.histogram2d(x[n], y[n], weights=key[n], bins=bins, range = ((-30, 30), (-30,30)))
        #histform, xbins, ybins = np.histogram2d(x_s[n], y_s[n], weights=tempform[n], bins=400, range = ((-30, 30), (-30,30)))
        im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='coolwarm', vmin = 2.1, vmax = 6)
        ax.set_xlim(-x_y_lim, x_y_lim)
        ax.set_ylim(-x_y_lim, x_y_lim)
        ax.text(0.5, 0.88, titlelist[n], horizontalalignment='center', transform=ax.transAxes)
        ax.set_xticklabels([])
        if (n == len(model)-1):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(T) [K]', size=12)
        if (n == 0):
            ax.set_ylabel('y [kpc]', fontsize = 12)

        else:
            ax.set_yticklabels([])

    # side-on
    else:
        ax = fig.add_subplot(gs0[n])
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(90)
        hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=bins, range = ((-30, 30), (-30,30)))
        im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='seismic', transform = rot+base, vmin = 2.1, vmax = 6)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlim(-x_y_lim, x_y_lim)
        ax.set_ylim(-z_lim, z_lim)
        ax.set_xlabel('x [kpc]', fontsize = 12)        
        if (n == 2*len(model)-1):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical')  
            
        if (n == len(model)):   
            ax.set_ylabel('z [kpc]', fontsize = 12)
        
        else:
            ax.set_yticklabels([])


fig.suptitle('Gas temperature')
plt.savefig('temp_map_' + model_name + '.pdf', bbox_inches='tight')
plt.clf()

