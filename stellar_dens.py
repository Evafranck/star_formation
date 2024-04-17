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
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fancybox'] = True
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.shadow'] = False
plt.rcParams['legend.edgecolor'] = 'darkgray'
plt.rcParams['patch.linewidth'] = 0.5

#young = f.LowPass('age', '1.49 Gyr')
#young = f.LowPass('age', '100 Myr')
bins = 400

density = []
x = []
y = []
x_y_lim = 14.99
z_lim = x_y_lim * 0.266

model_name = 'federrath_padoan' # options: 'hopkins', 'semenov', 'federrath_padoan', '1e6'
if model_name == '1e6':
    #models with 10**6 resolution
    model = ['threshold_1e6_alpha008', 'padoan_1e6_alpha008', 'federrath_1e6_alpha008', 'federrath_cstar_cut_1e6', 'semenov_1e6_alpha008', 'semenov_cstar_cut_1e6', 'evans_1e6_alpha008']

elif model_name == 'hopkins':
    # models that are based on Hopkins et al. (2013)
    model = ['threshold', 'hopkins', 'hopkins_alpha', 'hopkins_alpha008', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']

elif model_name == 'semenov':
    # models that are based on Semenov et al. (2016)
    model = ['threshold_alpha008', 'semenov_alpha008', 'semenov_cstar_cut', 'semenov_alpha008_tcr', 'semenov_alpha008_converging_flow', 'semenov_alpha008_tcool_cut', 'semenov_alpha008_tcool_cut_converging_flow']

elif model_name == 'federrath_padoan':
    # models that are based on Federrath & Klessen (2012) and Padoan et al. (2012)
    model = ['threshold_alpha008', 'federrath_alpha008', 'federrath_cstar_cut', 'padoan_alpha008']

titlelist = model

def load_sim_faceon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    #s = s_all.s[young]
    s = s_all.s
    s['n'] = s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s['n'])
    x.append(s['x'])
    y.append(s['y'])
    #print(mod, s.g['rho'].min(), s.g['rho'].max(), np.highian(s.g['rho']))
    print(mod, 'dens_min = ', s['n'].min(), 'dens_max = ', s['n'].max(), 'dens_mean = ', np.mean(s['n']))

def load_sim_sideon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    #s = s_all.s[young]
    s = s_all.s
    s['n'] = s['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    density.append(s['n'])
    x.append(s['x'])
    y.append(s['y'])
    #print(mod, s.g['rho'].min(), s.g['rho'].max(), np.median(s.g['rho']))
    print(mod, 'dens_min = ', s['n'].min(), 'dens_max = ', s['n'].max(), 'dens_mean = ', np.mean(s['n']))
    
for m in model:
    load_sim_faceon(m)
for m in model:    
    load_sim_sideon(m)
    
fig = plt.figure(figsize = (len(model)*2.33, 2.95))
width = np.ones(len(model))
width[-1] = 1.077
gs0 = gd.GridSpec(2, len(model), figure=fig, height_ratios = [1, 0.27], width_ratios = width)
gs0.update(hspace=0.00, wspace=0.00)

for n in range(len(model)):
    ax = fig.add_subplot(gs0[n])
    #print(s2.s['n'].max())
    #print(s2.s['n'].min())
    #print(s2.s['x'].max())
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], density[n], statistic='mean', bins=bins, range = ((-30, 30), (-30,30)))
    #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=density[n], bins=600, range = ((-50, 50), (-50,50)))
    im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r', vmin = -1.9, vmax = 2)
    #im = ax.imshow(np.log10(hist), extent=(-50,50,-50,50), cmap='CMRmap_r', vmin = 0.1, vmax = 5)
    if n == len(model)-1:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(n) [particles $cm^{-3}$]', size=10)
    
    if n == 0:    
        ax.set_ylabel('y [kpc]', fontsize = 10)
    else:
        ax.set_yticklabels([])
    
    ax.text(0.5, 0.88, titlelist[n], fontsize = 8, horizontalalignment='center', transform=ax.transAxes)
    ax.set_xlabel('x [kpc]', fontsize = 14)
    ax.set_xlim(-x_y_lim, x_y_lim)
    ax.set_ylim(-x_y_lim, x_y_lim)
    ax.set_aspect(1./ax.get_data_ratio())

for n in range(len(model), len(model)*2):
    ax = fig.add_subplot(gs0[n])
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(90)
    #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=density[n], bins=600, range = ((-30, 30), (-30,30)))
    #im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r', transform = rot+base, vmin = 0, vmax = 5)
    hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], density[n], statistic='mean', bins=bins, range = ((-30, 30), (-30,30)))
    im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r', transform = rot+base, vmin = -1.9, vmax = 2)
    
    if n == len(model)*2-1:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size = '5%', pad = 0.05)
        fig.colorbar(im, cax = cax, orientation='vertical')
    if n == len(model):
        ax.set_ylabel('z [kpc]', fontsize = 10)
    else:
        ax.set_yticklabels([])

    #ax.set_title(titlelist[n], fontsize = 11)
    ax.set_xlabel('x [kpc]', fontsize = 10)
    ax.set_xlim(-x_y_lim, x_y_lim)
    ax.set_ylim(-z_lim, z_lim)

plt.savefig('stellar_dens_' + model_name + '.pdf', bbox_inches='tight')
plt.clf()

