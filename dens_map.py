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

model_name = 'hopkins' # options: 'hopkins', 'semenov', 'federrath_padoan', '1e6'
key = []
x = []
y = []
range_tuple = ((-40, 40), (-40,40))
x_y_lim = 19.99
z_lim = x_y_lim * 0.266

def load_sim_faceon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s_disk = s_all[disk]
    cold = f.LowPass('temp', '30000 K') # nur kaltes gas
    s = s_disk.g[cold]
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    key.append(s.g['n'])
    x.append(s.g['x'])
    y.append(s.g['y'])
    print(mod, 'number of particles = ', len(s['x']), 'dens_min = ', s.g['n'].min(), 'dens_max = ', s.g['n'].max(), 'dens_mean = ', np.mean(s.g['n']))

def load_sim_sideon(mod):
    s_all = pynbody.load('../'+mod+'/halo.00128')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s_disk = s_all[disk]
    cold = f.LowPass('temp', '30000 K') # nur kaltes gas
    s = s_disk.g[cold]
    s.g['n'] = s.g['rho'].in_units('kg cm^-3')/(1.673*10**(-27))
    key.append(s.g['n'])
    x.append(s.g['x'])
    y.append(s.g['y'])

if model_name == '1e6':
    #models with 10**6 resolution
    model = ['threshold_1e6_alpha008', 'padoan_1e6_alpha008', 'federrath_1e6_alpha008', 'federrath_cstar_cut_1e6', 'semenov_1e6_alpha008', 'semenov_cstar_cut_1e6', 'evans_1e6_alpha008']
    bins = 400
elif model_name == 'hopkins':
    # models that are based on Hopkins et al. (2013)
    model = ['threshold', 'hopkins', 'hopkins_alpha', 'hopkins_alpha008', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
    bins = 230
elif model_name == 'semenov':
    # models that are based on Semenov et al. (2016)
    model = ['threshold_alpha008', 'semenov_alpha008', 'semenov_cstar_cut', 'semenov_alpha008_tcr', 'semenov_alpha008_converging_flow', 'semenov_alpha008_tcool_cut', 'semenov_alpha008_tcool_cut_converging_flow']
    bins = 230
elif model_name == 'federrath_padoan':
    # models that are based on Federrath & Klessen (2012) and Padoan et al. (2012)
    model = ['threshold_alpha008', 'federrath_alpha008', 'federrath_cstar_cut', 'padoan_alpha008']
    bins = 260
    
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
        #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=400, range = ((-30, 30), (-30,30)))
        hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], key[n], statistic='mean', bins=bins, range = range_tuple)
        #im = ax.imshow(np.log10((hist/surface_faceon)/(4*units.kpc).in_units('cm')), extent=(-30,30,-30,30), cmap='CMRmap_r')#, vmin = -2, vmax = 2)
        im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r') #, vmin = -1.9, vmax = 2)
        ax.set_xlim(-x_y_lim, x_y_lim)
        ax.set_ylim(-x_y_lim, x_y_lim)
        ax.text(0.5, 0.88, titlelist[n], fontsize = 8, horizontalalignment='center', transform=ax.transAxes)
        ax.set_xticklabels([])
        
        if (n == len(model)-1):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical').set_label(label = r'log(n) [particles $\mathrm{cm}^{-3}$]', size=12)
        
        if (n == 0):
            ax.set_ylabel('y [kpc]', fontsize = 10)
        else:
            ax.set_yticklabels([])

    # side-on
    else: 
        ax = fig.add_subplot(gs0[n])
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(90)
        hist, xbin, ybin, binnum = scipy.stats.binned_statistic_2d(x[n], y[n], key[n], statistic='mean', bins=bins, range = range_tuple)
        #hist, xbin, ybin = np.histogram2d(x[n], y[n],weights=key[n], bins=400, range = ((-30, 30), (-30,30)))
        #im = ax.imshow(np.log10((hist/surface_sideon)/(360*units.kpc).in_units('cm')), extent=(-30,30,-30,30), cmap='CMRmap_r', transform = rot+base)#, vmin = -2, vmax = 2)
        im = ax.imshow(np.log10(hist), extent=(-30,30,-30,30), cmap='CMRmap_r', transform = rot+base) #, vmin = -2, vmax = 2)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlim(-x_y_lim, x_y_lim)
        ax.set_ylim(-z_lim, z_lim)
        ax.set_xlabel('x [kpc]', fontsize = 10)        

        if (n == 2*len(model)-1):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0.05)
            fig.colorbar(im, cax = cax, orientation='vertical')  
            
        if (n == len(model)):
            ax.set_ylabel('z [kpc]', fontsize = 10)
        else:
            ax.set_yticklabels([])


fig.suptitle('Gas density')
plt.savefig('dens_map_' + model_name +'.pdf', bbox_inches='tight')
plt.clf()

