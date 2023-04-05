# function definitions to fit surface density

import numpy as np
import matplotlib.pylab as plt
import pynbody as pb
import os, pickle
import pynbody.filt as filt
import matplotlib.gridspec as gridspec
import scipy

from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import astropy 

# Define model
from astropy.modeling import Fittable1DModel, Parameter
from scipy.optimize import curve_fit
plt.switch_backend('agg') 
#%matplotlib inline
#%config InlineBackend.figure_format='retina'

import seaborn as sns


import pynbody
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
from pynbody import units
import pynbody.filt as f
from pynbody.analysis import profile
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
import astropy
import os, pickle
import seaborn as sns
# Define model
from astropy.modeling import Fittable1DModel, Parameter

#####################################################################
path = ['low_master_iso', 'low_padoan_iso', 'low_semenov_iso', 'low_evans_iso']
path2 = ['med_master_iso', 'med_padoan_iso', 'med_semenov_iso', 'med_evans_iso']
simulation = ['Threshold-based model, Low Resolution', 'Padoan et al. (2012), Low Resolution', 'Semenov et al. (2016), Low Resolution', 'Evans et al. (2022), Low Resolution', 'Threshold-based model, Medium Resolution', 'Padoan et al. (2012), Medium Resolution', 'Semenov et al. (2016), Medium Resolution', 'Evans et al. (2022), Medium Resolution']
#####################################################################

def exponential(x, amplitude, r_d):
    return np.log10(amplitude*np.exp(-x / r_d))


fig = plt.figure(figsize=(20,10))
gs = gridspec.GridSpec(2, 4)
gs.update(hspace=0.00, wspace=0.00)

axes = []
for i, cell in enumerate(gs):
    axes.append(plt.subplot(cell))
    if i == 1:
        axes[-1].set_yticklabels([])
    if i == 2:
        axes[-1].set_yticklabels([])
    if i == 3:
        axes[-1].set_yticklabels([])
    if i == 5:
        axes[-1].set_yticklabels([])
    if i == 6:
        axes[-1].set_yticklabels([])
    if i == 7:
        axes[-1].set_yticklabels([])
    if i == 0:
        axes[-1].set_xticklabels([])

for k, sim in enumerate(path):

    if not os.path.isfile(sim+'_surf_den.dat'):
        s = pb.load(sim+'/low.01000')
        s.physical_units()
        new = filt.LowPass('age', s.properties['time'].in_units('Gyr'))
        h1 = s.s[new]
        #h = s.halos()
        #h1 = h[1]
        #pb.analysis.angmom.faceon(h1)

        ''''''''''
        mass_star = (h1.s['mass'].in_units('Msol')).view(np.ndarray)
        rxy_star = np.sqrt((h1.s['x']).view(np.ndarray)**2 + (h1.s['y']).view(np.ndarray)**2)
        srt = np.argsort(rxy_star)
        mass_star = mass_star[srt]
        rxy_star = rxy_star[srt]
    
        rgal = 30 #kpc  #2 * r_90[k]
        eps = h1['eps'].min()
        nbins = int(rgal/eps)
        
        r = rxy_star[rxy_star<=rgal]
        m = mass_star[rxy_star<=rgal]

        every_n = int(len(r)/nbins)
        bins_edges = r[::every_n]
        bin_mass , bin_edges, bin_number = scipy.stats.binned_statistic(r, m, statistic='sum', bins=bins_edges)
        bin_massr , bin_edges, bin_number = scipy.stats.binned_statistic(r, m*r, statistic='sum', bins=bins_edges)
        R_bins = bin_massr/bin_mass # then bins' ce:wqnters are computed as mass weigthed radii
        A_bins = np.array([np.pi*(bin_edges[i]**2-bin_edges[i-1]**2) for i in range(1,len(bin_edges))]) # the bins' areas are in kpc**2

        log_surface_mass_density = np.log10(bin_mass/A_bins)-6.0 # units log(Msol/pc**2) 
        surface_mass_density = bin_mass/(A_bins) # units Msol/kpc**2 
        #fitting sersic plus exponential
        g = fit_g(s_init, R_bins, log_surface_mass_density)
        #print(g)
        #surf = fit_g(s_init, R_bins, surface_mass_density)
        #print surf
        exp1 = exponential(g.amplitude, g.r_d)
        '''''''''''
        ps_new = profile.VerticalProfile(h1, '1.35 kpc', '10 kpc', '3 kpc', ndim = 2, nbins = 20)
        
        R_bins = ps_new['rbins'].view(np.ndarray)[7:-1]
        ps_surf_dens = ps_new['density'].view(np.ndarray)[7:-1]
        print(ps_surf_dens)
        popt, pcov = curve_fit(exponential, R_bins, np.log10(ps_surf_dens), bounds=(0.0001, [1e15, 5]))
        #log_ps_surf_dens = np.log10(ps_surf_dens)

        #f = open(sim+'_surf_den.dat','wb')
        #pickle.dump({'rbins':R_bins,'dens':log_surface_mass_density,'fit':g,'exp':exp1},f)
        #f.close()i
        #exp1 = models.Exponential1D(popt[0],popt[1])

    else:
        #load from pickle

        data = pickle.load(open(sim+'_surf_den.dat','rb'))
        R_bins = data['rbins']
        log_surface_mass_density = data['dens']
        g = data['fit']
        exp1 = data['exp']

    #plot the data
    axes[k].scatter(ps_new['rbins'].view(np.ndarray), np.log10(ps_new['density'].view(np.ndarray)),marker='o', s=1.5, color='k',zorder=10, label = 'simulation data')
    
    # plot exponential
    a = popt[0]/(10**9)
    sigma_a = np.sqrt(pcov[0][0])/(10**9)
    axes[k].plot(R_bins, exponential(R_bins, *popt),color='orange',lw=3,zorder=1, label = r'fit:'+"\n"+r"$R_{\rm h}=($%.3f"%popt[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov[1][1])+"\n"+r"$\Sigma_{\rm h}=($%.2f"%a+r'$\pm$%.2f)$ \cdot 10^9$ M$_{\rm\odot}$pc$^{\rm -2}$'%sigma_a)

    axes[k].set_xlim(-0.1, 4.4)
    axes[k].set_ylim(4.5, 10.4)
    if k == 0: axes[k].set_ylabel(r"log($\rm\Sigma$) [M$_{\rm\odot}$pc$^{\rm -2}$]", fontsize = 16)
    #axes[k].plot([r_50[k],r_50[k]],[0,4],ls='dotted',color='gray',lw=2)
    #axes[k].plot([r_90[k],r_90[k]],[0,4],ls='dashed',color='gray',lw=2)
    axes[k].legend(fontsize = 12, loc='upper center', bbox_to_anchor=(0.5, 0.3))
     
    # label panels
    axes[k].text(0.5,0.9,simulation[k],fontsize=12,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[k].transAxes)
    #axes[1].text(0.5,0.9,'g7.66e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[1].transAxes)
    #axes[2].text(0.5,0.9,'g8.26e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[2].transAxes)
    #axes[3].text(0.5,0.9,'g1.12e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[3].transAxes)
    #axes[4].text(0.5,0.9,'g2.79e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[4].transAxes)

for k, sim in enumerate(path2):

    if not os.path.isfile(sim+'_surf_den.dat'):
        s = pb.load(sim+'/med.01000')
        s.physical_units()
        new = filt.LowPass('age', s.properties['time'].in_units('Gyr'))
        h1 = s.s[new]
        #h = s.halos()
        #h1 = h[1]
        #pb.analysis.angmom.faceon(h1)

        ''''''''''
        mass_star = (h1.s['mass'].in_units('Msol')).view(np.ndarray)
        rxy_star = np.sqrt((h1.s['x']).view(np.ndarray)**2 + (h1.s['y']).view(np.ndarray)**2)
        srt = np.argsort(rxy_star)
        mass_star = mass_star[srt]
        rxy_star = rxy_star[srt]
    
        rgal = 30 #kpc  #2 * r_90[k]
        eps = h1['eps'].min()
        nbins = int(rgal/eps)
        
        r = rxy_star[rxy_star<=rgal]
        m = mass_star[rxy_star<=rgal]

        every_n = int(len(r)/nbins)
        bins_edges = r[::every_n]
        bin_mass , bin_edges, bin_number = scipy.stats.binned_statistic(r, m, statistic='sum', bins=bins_edges)
        bin_massr , bin_edges, bin_number = scipy.stats.binned_statistic(r, m*r, statistic='sum', bins=bins_edges)
        R_bins = bin_massr/bin_mass # then bins' ce:wqnters are computed as mass weigthed radii
        A_bins = np.array([np.pi*(bin_edges[i]**2-bin_edges[i-1]**2) for i in range(1,len(bin_edges))]) # the bins' areas are in kpc**2

        log_surface_mass_density = np.log10(bin_mass/A_bins)-6.0 # units log(Msol/pc**2) 
        surface_mass_density = bin_mass/(A_bins) # units Msol/kpc**2 
        #fitting sersic plus exponential
        g = fit_g(s_init, R_bins, log_surface_mass_density)
        #print(g)
        #surf = fit_g(s_init, R_bins, surface_mass_density)
        #print surf
        exp1 = exponential(g.amplitude, g.r_d)
        '''''''''''
        ps_new = profile.VerticalProfile(h1, '1.35 kpc', '10 kpc', '3 kpc', ndim = 2, nbins = 20)
        
        if sim == 'med_master_iso':
            R_bins1 = ps_new['rbins'].view(np.ndarray)[:7]
            ps_surf_dens1 = ps_new['density'].view(np.ndarray)[:7]
            R_bins2 = ps_new['rbins'].view(np.ndarray)[7:-1]
            ps_surf_dens2 = ps_new['density'].view(np.ndarray)[7:-1]
            print(ps_surf_dens2)
            popt1, pcov1 = curve_fit(exponential, R_bins1, np.log10(ps_surf_dens1), bounds=(0.0001, [1e15, 5]))
            popt2, pcov2 = curve_fit(exponential, R_bins2, np.log10(ps_surf_dens2), bounds=(0.0001, [1e15, 5])) 

        if sim == 'med_semenov_iso':
            R_bins1 = ps_new['rbins'].view(np.ndarray)[:7]
            ps_surf_dens1 = ps_new['density'].view(np.ndarray)[:7]
            R_bins2 = ps_new['rbins'].view(np.ndarray)[7:]
            ps_surf_dens2 = ps_new['density'].view(np.ndarray)[7:]
            print(ps_surf_dens1)
            print(ps_surf_dens2)
            popt1, pcov1 = curve_fit(exponential, R_bins1, np.log10(ps_surf_dens1), bounds=(0.0001, [1e15, 5]))
            popt2, pcov2 = curve_fit(exponential, R_bins2, np.log10(ps_surf_dens2), bounds=(0.0001, [1e15, 5]))

        if sim == 'med_padoan_iso':
            R_bins = ps_new['rbins'].view(np.ndarray)[:-1]
            ps_surf_dens = ps_new['density'].view(np.ndarray)[:-1]
            print(ps_surf_dens)
            popt, pcov = curve_fit(exponential, R_bins, np.log10(ps_surf_dens), bounds=(0.0001, [1e15, 5]))
        
        if sim == 'med_evans_iso':
            R_bins = ps_new['rbins'].view(np.ndarray)[:-1]
            ps_surf_dens = ps_new['density'].view(np.ndarray)[:-1]
            print(ps_surf_dens)
            popt, pcov = curve_fit(exponential, R_bins, np.log10(ps_surf_dens), bounds=(0.0001, [1e15, 5]))
        #log_ps_surf_dens = np.log10(ps_surf_dens)

        #f = open(sim+'_surf_den.dat','wb')
        #pickle.dump({'rbins':R_bins,'dens':log_surface_mass_density,'fit':g,'exp':exp1},f)
        #f.close()i
        #exp1 = models.Exponential1D(popt[0],popt[1])

    else:
        #load from pickle

        data = pickle.load(open(sim+'_surf_den.dat','rb'))
        R_bins = data['rbins']
        log_surface_mass_density = data['dens']
        g = data['fit']
        exp1 = data['exp']

    #plot the data
    axes[k+4].scatter(ps_new['rbins'].view(np.ndarray), np.log10(ps_new['density'].view(np.ndarray)),marker='o', s=1.5, color='k',zorder=10, label = 'simulation data')
    
    #axes[k].plot(R_bins,g(R_bins),color='blue',lw=3,zorder=0, label = r'fit for the disk:'+"\n"+r"$\Sigma_{\rm 0}=$%.2f M$_{\rm\odot}$pc$^{\rm -2}$"%g.amplitude.value+"\n"+r"$R_{\rm d}=$%.2f kpc"%g.r_d.value) 
    
    # plot exponential
    if sim == 'med_master_iso':
        a1 = popt1[0]/(10**9)
        sigma_a1 = np.sqrt(pcov1[0][0])/(10**9)
        a2 = popt[0]/(10**9)
        sigma_a2 = np.sqrt(pcov2[0][0])/(10**9)
        axes[k+4].plot(R_bins1, exponential(R_bins1, *popt1),color='orange',lw=3,zorder=1, label = r'fit:'+"\n"+r"$R_{\rm h, 1}=($%.3f"%popt1[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov1[1][1])+"\n"+r"$\Sigma_{\rm h, 1}=($%.2f"%a1+r'$\pm$%.2f)$ \cdot 10^9$ M$_{\rm\odot}$pc$^{\rm -2}$'%sigma_a1)
        axes[k+4].plot(R_bins2, exponential(R_bins2, *popt2),color='red',lw=3,zorder=1, label = r'fit:'+"\n"+r"$R_{\rm h, 2}=($%.3f"%popt2[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov2[1][1])+"\n"+r"$\Sigma_{\rm h, 2}=($%.2f"%a2+r'$\pm$%.2f)$ \cdot 10^9$ M$_{\rm\odot}$pc$^{\rm -2}$'%sigma_a2)
        axes[k+4].legend(fontsize = 10, loc='upper center', bbox_to_anchor=(0.64, 0.93))

    if sim == 'med_semenov_iso':
        a1 = popt1[0]/(10**9)
        sigma_a1 = np.sqrt(pcov1[0][0])/(10**9)
        a2 = popt2[0]/(10**9)
        sigma_a2 = np.sqrt(pcov2[0][0])/(10**9)
        axes[k+4].plot(R_bins1, exponential(R_bins1, *popt1),color='orange',lw=3,zorder=1, label = r'fit:'+"\n"+r"$R_{\rm h, 1}=($%.3f"%popt1[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov1[1][1])+"\n"+r"$\Sigma_{\rm h, 1}=($%.2f"%a1+r'$\pm$%.2f)$ \cdot 10^9$ M$_{\rm\odot}$pc$^{\rm -2}$'%sigma_a1)
        axes[k+4].plot(R_bins2, exponential(R_bins2, *popt2),color='red',lw=3,zorder=1, label = r'fit:'+"\n"+r"$R_{\rm h, 2}=($%.3f"%popt2[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov2[1][1])+"\n"+r"$\Sigma_{\rm h, 2}=($%.2f"%a2+r'$\pm$%.2f)$ \cdot 10^9$ M$_{\rm\odot}$pc$^{\rm -2}$'%sigma_a2)
        axes[k+4].legend(fontsize = 10, loc='upper center', bbox_to_anchor=(0.65, 0.93))

    if sim == 'med_padoan_iso':
        a = popt[0]/(10**9)
        sigma_a = np.sqrt(pcov[0][0])/(10**9)
        axes[k+4].plot(R_bins, exponential(R_bins, *popt),color='orange',lw=3,zorder=1, label = r'fit:'+"\n"+r"$R_{\rm h}=($%.3f"%popt[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov[1][1])+"\n"+r"$\Sigma_{\rm h}=($%.2f"%a+r'$\pm$%.2f)$ \cdot 10^9$ M$_{\rm\odot}$pc$^{\rm -2}$'%sigma_a)
        axes[k+4].legend(fontsize = 11, loc='upper center', bbox_to_anchor=(0.6, 0.9))

    if sim == 'med_evans_iso':
        a = popt[0]/(10**9)
        sigma_a = np.sqrt(pcov[0][0])/(10**9)
        axes[k+4].plot(R_bins, exponential(R_bins, *popt),color='orange',lw=3,zorder=1, label = r'fit:'+"\n"+r"$R_{\rm h}=($%.3f"%popt[1]+r'$\pm$%.3f) kpc'%np.sqrt(pcov[1][1])+"\n"+r"$\Sigma_{\rm h}=($%.2f"%a+r'$\pm$%.2f)$ \cdot 10^9$ M$_{\rm\odot}$pc$^{\rm -2}$'%sigma_a)
        axes[k+4].legend(fontsize = 11, loc='upper center', bbox_to_anchor=(0.6, 0.9))

    axes[k+4].set_xlim(-0.1, 4.4)
    axes[k+4].set_ylim(4.5, 10.4)
    axes[k+4].set_xlabel(r"z [kpc]", fontsize = 16)
    if k == 0: axes[k+4].set_ylabel(r"log($\rm\Sigma$) [M$_{\rm\odot}$pc$^{\rm -2}$]", fontsize = 16)
    #axes[k+4].plot([r_50[k+4],r_50[k+4]],[0,4],ls='dotted',color='gray',lw=2)
    #axes[k+4].plot([r_90[k+4],r_90[k+4]],[0,4],ls='dashed',color='gray',lw=2)
     
    # label panels
    axes[k+4].text(0.5,0.95,simulation[k+4],fontsize=12,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[k+4].transAxes)
    #axes[1].text(0.5,0.9,'g7.66e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[1].transAxes)
    #axes[2].text(0.5,0.9,'g8.26e11',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[2].transAxes)
    #axes[3].text(0.5,0.9,'g1.12e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[3].transAxes)
    #axes[4].text(0.5,0.9,'g2.79e12',fontsize=35,color='black',horizontalalignment='center',verticalalignment='center',transform=axes[4].transAxes)


plt.savefig('surface_density_fits_vert_both.pdf', bbox_inches='tight')
