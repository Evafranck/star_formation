import pynbody
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gd
from matplotlib import transforms
from pynbody import units
from pynbody.analysis import profile
import pynbody.filt as f
from pynbody import array
import glob
import scipy.stats

plt.rc('axes', linewidth=1)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 1
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 1
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['axes.edgecolor'] = 'k'#'gray'
#plt.rcParams['axes.grid'] = True
#plt.rcParams['grid.color'] = 'lightgray'
#plt.rcParams['grid.linestyle'] = 'dashed' #dashes=(5, 1)
plt.rcParams['lines.dashed_pattern'] = 10, 3
plt.rcParams['grid.linewidth'] = 1.5
#plt.rcParams['axes.facecolor'] = 'whitesmoke'
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['legend.fancybox'] = True
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.shadow'] = False
plt.rcParams['legend.edgecolor'] = 'darkgray'
plt.rcParams['patch.linewidth'] = 1

rbins = []
zbins = []
rbins_gas = []
zbins_gas = []
star_surf_dens = []
gas_surf_dens = []
star_vert_surf_dens = []
gas_vert_surf_dens = []


#titlelist = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r'$\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)']
#simulations = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
simulations = ['threshold_alpha008','threshold_1e6_alpha008', 'semenov_alpha008', 'semenov_1e6_alpha008', 'semenov_cstar_cut', 'federrath_alpha008', 'federrath_1e6_alpha008', 'federrath_cstar_cut'] 
#simulations = ['threshold_alpha008', 'threshold_1e6_alpha008', 'hopkins_alpha_padoan_alpha008', 'hopkins_alpha008', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
titlelist = simulations

# Load a slice of the simulation snapshots faceon and sideon
def load_sim_faceon(mod):
    s_all = pynbody.load('../' + mod + '/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    #new = f.LowPass('age', s_all.properties['time'].in_units('Gyr'))
    new = f.LowPass('age', '1.49 Gyr')
    s_new = s_all.s[new]
    cold = f.LowPass('temp', '30000 K') # nur kaltes gas
    s_gas = s_all.g[cold]
    ps = pynbody.analysis.profile.Profile(s_new, nbins=20, rmin='0 kpc', rmax = '25 kpc')
    ps_vert = profile.VerticalProfile(s_new, '0.1 kpc', '25 kpc', '5 kpc', ndim = 2, nbins = 20)
    ps_gas = pynbody.analysis.profile.Profile(s_gas, nbins=20, rmin='0 kpc', rmax = '25 kpc')
    ps_vert_gas = profile.VerticalProfile(s_gas, '0.1 kpc', '25 kpc', '5 kpc', ndim = 2, nbins = 20)
    rbins.append(ps['rbins'].in_units('kpc'))
    zbins.append(ps_vert['rbins'].in_units('kpc'))
    star_surf_dens.append(np.log10(ps['density'].in_units('Msol pc^-2')))
    star_vert_surf_dens.append(np.log10(ps_vert['density'].in_units('Msol pc^-2')))
    rbins_gas.append(ps_gas['rbins'].in_units('kpc'))
    zbins_gas.append(ps_vert_gas['rbins'].in_units('kpc'))
    gas_surf_dens.append(np.log10(ps_gas['density'].in_units('Msol pc^-2')))
    gas_vert_surf_dens.append(np.log10(ps_vert_gas['density'].in_units('Msol pc^-2')))

    
for m in simulations:
    load_sim_faceon(m)
    

fig = plt.figure(figsize = (23,5))
gs0 = gd.GridSpec(1, 4, figure=fig, width_ratios = [1, 1, 1, 1])

for n in range(4):
    ax = fig.add_subplot(gs0[n])  
    if (n==0):
        for i in range(8):
            plt.plot(rbins[i], star_surf_dens[i],  label = titlelist[i], lw=1) #, c = colorlist[i], ls = linestyle[i])
        ax.set_xlabel('R [kpc]', fontsize = 14)
        #ax.set_xlim(0, 19)
        #ax.set_ylim(-4, 3)
        ax.legend(fontsize = 10, loc = 'upper right')
        ax.set_ylabel(r'log($\Sigma_{\star}$) [M$_{\odot}$ kpc$^{-2}$]', fontsize = 14)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_title('a) Radial surface density profile of stars', fontsize = 14)
        
    if (n==1):
        for i in range(8):
            plt.plot(rbins_gas[i], gas_surf_dens[i], lw = 1) #, c = colorlist[i], ls = linestyle[i])
        ax.set_xlabel('R [kpc]', fontsize = 14)
        #ax.set_xlim(0, 4)
        #ax.set_ylim(-3, 3)
        ax.set_ylabel(r'log($\Sigma_{gas}$) [M$_{\odot}$ kpc$^{-2}$]', fontsize = 14)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_title('b) Radial surface density profile of gas', fontsize = 14)
        
    if (n==2):
        for i in range(8):
            plt.plot(zbins_gas[i], gas_vert_surf_dens[i], lw = 1) #, c = colorlist[i], ls = linestyle[i])
        ax.set_xlabel('z [kpc]', fontsize = 14)
        #ax.set_xlim(0, 4)
        #ax.set_ylim(-3, 3)
        ax.set_ylabel(r'log($\Sigma_{gas}$) [M$_{\odot}$ kpc$^{-2}$]', fontsize = 14)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_title('c) Vertical surface density profile of gas', fontsize = 14)
    
    if (n==3):
        for i in range(6):
            plt.plot(zbins_gas[i], star_vert_surf_dens[i], lw = 1) #, c = colorlist[i], ls = linestyle[i])
        ax.set_xlabel('z [kpc]', fontsize = 14)
        #ax.set_xlim(0, 4)
        #ax.set_ylim(-3, 3)
        ax.set_ylabel(r'log($\Sigma_{gas}$) [M$_{\odot}$ kpc$^{-2}$]', fontsize = 14)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_title('d) Vertical surface density profile of stars', fontsize = 14)
        
plt.savefig('surf_dens2.pdf', bbox_inches = 'tight')



