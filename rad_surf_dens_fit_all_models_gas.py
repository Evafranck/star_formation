# function definitions to fit surface density

import numpy as np
import matplotlib.pylab as plt
import pynbody as pb
import os, pickle

import matplotlib.gridspec as gridspec
import scipy
import pynbody.filt as filt

from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import astropy 

# Define model
from astropy.modeling import Fittable1DModel, Parameter

from scipy.optimize import curve_fit

#plt.switch_backend('agg') 
#%matplotlib inline
#%config InlineBackend.figure_format='retina'

import seaborn as sns
plt.rcParams['figure.figsize'] = (15, 10)

sns.set_style('ticks')
#sns.set_style('darkgrid')
sns.set_context("talk",font_scale=1,rc={"lines.linewidth": 2,"axes.linewidth": 2})

plt.rc('axes', linewidth=2)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 3
plt.rcParams['ytick.minor.size'] = 4
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
plt.rcParams['patch.linewidth'] = 2


def surface_density(r, amplitude1, r_eff, n, amplitude2, r_d):
    from scipy.special import gammaincinv
    return (np.log10(amplitude1 * np.exp(-gammaincinv(2 * n, 0.5) * ((r / r_eff) ** (1 / n) - 1)) + amplitude2 * np.exp(- r / r_d)))

#####################################################################################################

path = ['low_master_iso', 'low_padoan_iso', 'low_semenov_iso', 'low_evans_iso']
path2 = ['med_master_iso', 'med_padoan_iso', 'med_semenov_iso', 'med_evans_iso']
simulation = ['Threshold-based model'+'\n'+'Low Resolution', 'Padoan et al. (2012)'+'\n'+'Low Resolution','Semenov et al. (2016)'+'\n'+'Low Resolution', 'Evans et al. (2022)'+'\n'+'Low Resolution', 'Threshold-based model'+'\n'+'Medium Resolution', 'Padoan et al. (2012)'+'\n'+'Medium Resolution','Semenov et al. (2016)'+'\n'+'Medium Resolution', 'Evans et al. (2022)'+'\n'+'Medium Resolution']

######################################################################################################

fig = plt.figure(figsize=(20,10))
gs = gridspec.GridSpec(2, 4)
gs.update(hspace=0.00, wspace=0.00)

#fig.suptitle('Surface density fit', fontsize = 16)
axes = []
for i, cell in enumerate(gs):
    axes.append(plt.subplot(cell))
    if i == 1:
        axes[-1].set_yticklabels([])
    if i == 5:
        axes[-1].set_yticklabels([])
    if i == 6:
        axes[-1].set_yticklabels([])
    if i == 3:
        axes[-1].set_yticklabels([])
    if i == 2:
        axes[-1].set_yticklabels([])
    if i == 7:
        axes[-1].set_yticklabels([])

for k, sim in enumerate(path):
    sim = '/mnt/storage/evafranck_data/simulations/agora/'+sim
    if not os.path.isfile(sim+'_surf_den_test.dat'):
        h1 = pb.load(sim + '/low.01000')
        h1.physical_units()
        #cold = filt.LowPass('temp', '30000 K')
        #h1 = s.g[cold]
        #h = s.halos()
        #h1 = h[1]
        #h1.physical_units()
        #pb.analysis.angmom.faceon(h1)

        mass_star = (h1.g['mass'].in_units('Msol')).view(np.ndarray)
        rxy_star = np.sqrt((h1.g['x']).view(np.ndarray)**2 + (h1.g['y']).view(np.ndarray)**2)
        srt = np.argsort(rxy_star)
        mass_star = mass_star[srt]
        rxy_star = rxy_star[srt]
    
        rgal = 30 #kpc  #2 * r_90[k]
        eps = 0.1#h1['eps'].min() #this gives 50 bins
        nbins = int(rgal/eps)
        
        r = rxy_star[rxy_star<=rgal]
        m = mass_star[rxy_star<=rgal]

        every_n = int(len(r)/nbins)
        bins_edges = r[::every_n]
        bin_mass , bin_edges, bin_number = scipy.stats.binned_statistic(r, m, statistic='sum', bins=bins_edges)
        bin_massr , bin_edges, bin_number = scipy.stats.binned_statistic(r, m*r, statistic='sum', bins=bins_edges)
        R_bins = bin_massr/bin_mass # then bins' centers are computed as mass weigthed radii
        A_bins = np.array([np.pi*(bin_edges[i]**2-bin_edges[i-1]**2) for i in range(1,len(bin_edges))]) # the bins' areas are in kpc**2

        log_surface_mass_density = np.log10(bin_mass/A_bins)-6.0 # units log(Msol/pc**2) 
        surface_mass_density = bin_mass/(A_bins)/1000000. # units Msol/kpc**2 
        
        popt, pcov = curve_fit(surface_density, R_bins, log_surface_mass_density, bounds=(0, [1e8, 2., 3., 1e8, 5]))
        #print(pcov)
        sersic = models.Sersic1D(popt[0],popt[1],popt[2])#g.amplitude_0, g.r_eff_0, g.n_0) 
        exp1 = models.Exponential1D(popt[3],popt[4])#g.amplitude_1, g.tau_1)

        #f = open(sim+'_surf_den.dat','wb')
        #pickle.dump({'rbins':R_bins,'dens':surface_mass_density,'fit':g,'sersic':sersic,'exp':exp1},f)
        #f.close()

    else:
        #load from pickle

        data = pickle.load(open(sim+'_surf_den.dat','rb'))
        R_bins = data['rbins']
        log_surface_mass_density = data['dens']
        g = data['fit']
        sersic = data['sersic']
        exp1 = data['exp']
        #print(g.fit_info['param_cov'])
    #plot the data
    axes[k].scatter(R_bins,np.log10(surface_mass_density),marker='o', s=1.5, color='k',zorder=10, label = 'simulation data')
    print(popt[4]) 
    # plot sersic
    #axes[k].plot(R_bins,np.log10(sersic(R_bins)),color='r',lw=3,zorder=1,label=r"fit for the bulge:"+"\n"+r"$n=$%.2f"%popt[2]+r'$\pm$%.2f'%np.sqrt(pcov[2][2])+"\n"+r"$R_{\rm eff}=($%.3f"%popt[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov[1][1])+'\n'+r"$\Sigma_{\rm eff}=($%.2f"%popt[0] +r'$\pm$%.2f) M$_{\rm\odot}$pc$^{\rm -2}$'%np.sqrt(pcov[0][0]))

    # plot sersic + exponential
    axes[k].plot(R_bins, surface_density(R_bins, *popt),color='orange',lw=3,zorder=1, label = r'fit for the disk:'+"\n"+r"$R_{\rm d}=($%.3f"%popt[4]+r'$\pm$%.3f) kpc'%np.sqrt(pcov[4][4])+"\n"+r"$\Sigma_{\rm d}=($%.2f"%popt[3]+r'$\pm$%.2f) M$_{\rm\odot}$pc$^{\rm -2}$'%np.sqrt(pcov[3][3]))


    axes[k].set_ylim(-1,3.4)
    axes[k].set_xlim(-1.5,29.9)
    if k == 0: axes[k].set_ylabel(r"log($\rm\Sigma$) [M$_{\rm\odot}$pc$^{\rm -2}$]", fontsize = 16)
    #axes[k].plot([r_50[k],r_50[k]],[0,4],ls='dotted',color='gray',lw=2)
    #axes[k].plot([r_90[k],r_90[k]],[0,4],ls='dashed',color='gray',lw=2)
    axes[k].legend(fontsize=13, loc='upper center', bbox_to_anchor=(0.58, 0.85))
     
    # label panels
    axes[k].text(0.5,0.9,simulation[k],fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[k].transAxes)
    #axes[1].text(0.5,0.9,'g7.66e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[1].transAxes)
    #axes[2].text(0.5,0.9,'g8.26e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[2].transAxes)
    #axes[3].text(0.5,0.9,'g1.12e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[3].transAxes)
    #axes[4].text(0.5,0.9,'g2.79e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[4].transAxes)

for k, sim in enumerate(path2):
    sim = '/mnt/storage/evafranck_data/simulations/agora/'+sim
    if not os.path.isfile(sim+'_surf_den_test.dat'):
        h1 = pb.load(sim + '/med.01000')
        h1.physical_units()
        #cold = filt.LowPass('temp', '30000 K')
        #h1 = s.g[cold]
        #h = s.halos()
        #h1 = h[1]
        #h1.physical_units()
        #pb.analysis.angmom.faceon(h1)

        mass_star = (h1.g['mass'].in_units('Msol')).view(np.ndarray)
        rxy_star = np.sqrt((h1.g['x']).view(np.ndarray)**2 + (h1.g['y']).view(np.ndarray)**2)
        srt = np.argsort(rxy_star)
        mass_star = mass_star[srt]
        rxy_star = rxy_star[srt]
    
        rgal = 30 #kpc  #2 * r_90[k]
        eps = 0.1#h1['eps'].min() #this gives 50 bins
        nbins = int(rgal/eps)
        
        r = rxy_star[rxy_star<=rgal]
        m = mass_star[rxy_star<=rgal]

        every_n = int(len(r)/nbins)
        bins_edges = r[::every_n]
        bin_mass , bin_edges, bin_number = scipy.stats.binned_statistic(r, m, statistic='sum', bins=bins_edges)
        bin_massr , bin_edges, bin_number = scipy.stats.binned_statistic(r, m*r, statistic='sum', bins=bins_edges)
        R_bins = bin_massr/bin_mass # then bins' centers are computed as mass weigthed radii
        A_bins = np.array([np.pi*(bin_edges[i]**2-bin_edges[i-1]**2) for i in range(1,len(bin_edges))]) # the bins' areas are in kpc**2

        log_surface_mass_density = np.log10(bin_mass/A_bins)-6.0 # units log(Msol/pc**2) 
        surface_mass_density = bin_mass/(A_bins)/1000000. # units Msol/kpc**2 
        
        popt, pcov = curve_fit(surface_density, R_bins, log_surface_mass_density, bounds=(0, [1e2, 2., 3., 1e2, 5]))
        #print(pcov)
        sersic = models.Sersic1D(popt[0],popt[1],popt[2])#g.amplitude_0, g.r_eff_0, g.n_0) 
        exp1 = models.Exponential1D(popt[3],popt[4])#g.amplitude_1, g.tau_1)

        #f = open(sim+'_surf_den.dat','wb')
        #pickle.dump({'rbins':R_bins,'dens':surface_mass_density,'fit':g,'sersic':sersic,'exp':exp1},f)
        #f.close()

    else:
        #load from pickle

        data = pickle.load(open(sim+'_surf_den.dat','rb'))
        R_bins = data['rbins']
        log_surface_mass_density = data['dens']
        g = data['fit']
        sersic = data['sersic']
        exp1 = data['exp']
        #print(g.fit_info['param_cov'])
    #plot the data
    print(popt[4])
    axes[k+4].scatter(R_bins[5:],np.log10(surface_mass_density[5:]),marker='o', s=1.5, color='k',zorder=10, label = 'simulation data')
    
    # plot sersic
    #axes[k+4].plot(R_bins,np.log10(sersic(R_bins)),color='r',lw=3,zorder=1,label=r"fit for the bulge:"+"\n"+r"$n=$%.2f"%popt[2]+r'$\pm$%.2f'%np.sqrt(pcov[2][2])+"\n"+r"$R_{\rm eff}=($%.3f"%popt[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov[1][1])+'\n'+r"$\Sigma_{\rm eff}=($%.2f"%popt[0] +r'$\pm$%.2f) M$_{\rm\odot}$pc$^{\rm -2}$'%np.sqrt(pcov[0][0]))

    # plot sersic + exponential
    axes[k+4].plot(R_bins, surface_density(R_bins, *popt),color='orange',lw=3,zorder=1, label = r'fit for the disk:'+"\n"+r"$R_{\rm d}=($%.3f"%popt[4]+r'$\pm$%.3f) kpc'%np.sqrt(pcov[4][4])+"\n"+r"$\Sigma_{\rm d}=($%.2f"%popt[3]+r'$\pm$%.2f) M$_{\rm\odot}$pc$^{\rm -2}$'%np.sqrt(pcov[3][3]))


    axes[k+4].set_ylim(-1,3.4)
    axes[k+4].set_xlim(-1.5,29.9)
    axes[k+4].set_xlabel('R [kpc]')
    if k == 0: axes[k+4].set_ylabel(r"log($\rm\Sigma$) [M$_{\rm\odot}$pc$^{\rm -2}$]", fontsize = 16)
    #axes[k+4].plot([r_50[k+4],r_50[k+4]],[0,4],ls='dotted',color='gray',lw=2)
    #axes[k+4].plot([r_90[k+4],r_90[k+4]],[0,4],ls='dashed',color='gray',lw=2)
    axes[k+4].legend(fontsize=13, loc='upper center', bbox_to_anchor=(0.58, 0.85))
     
    # label panels
    axes[k+4].text(0.5,0.9,simulation[k+4],fontsize=18,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[k+4].transAxes)
    #axes[1].text(0.5,0.9,'g7.66e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[1].transAxes)
    #axes[2].text(0.5,0.9,'g8.26e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[2].transAxes)
    #axes[3].text(0.5,0.9,'g1.12e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[3].transAxes)
    #axes[4].text(0.5,0.9,'g2.79e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[4].transAxes)


plt.savefig('surface_density_fits_both_gas.pdf', bbox_inches='tight')
