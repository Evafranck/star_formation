# function definitions to fit surface density

import numpy as np
import matplotlib.pylab as plt
import pynbody as pb
import os, pickle
import glob
from pathlib import Path
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

plt.rc('axes', linewidth=1)
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
plt.rcParams['patch.linewidth'] = 1


def surface_density(r, amplitude1, r_eff, n, amplitude2, r_d):
    from scipy.special import gammaincinv
    return (np.log10(amplitude1 * np.exp(-gammaincinv(2 * n, 0.5) * ((r / r_eff) ** (1 / n) - 1)) + amplitude2 * np.exp(- r / r_d)))

#####################################################################################################
pathlist = ['../high_master_iso', '../padoan', '../high_evans_iso', '../high_federrath_iso']
simulation = ['Threshold-based model', 'Padoan et al. (2012)', 'Evans et al. (2022)', 'Federrath & Klessen (2012)' + '\n' + 'without temperature cut']
######################################################################################################

fig = plt.figure(figsize=(10,10))
gs0 = gridspec.GridSpec(1, 1)
gs0.update(hspace=0.00, wspace=0.00)


label = ['Threshold-based model', 'Padoan et al. (2012)', 'Evans et al. (2022)', 'Federrath & Klessen (2012)']
label2 = ['bulge threshold-based model', 'bulge Semenov et al. (2016)', 'bulge Evans et al. (2022)', 'Federrath & Klessen (2012)']
color = ['blue', 'orange', 'green', 'red']
color2 = ['blue', 'orange', 'green', 'red']
for n in range(1):
    ax = fig.add_subplot(gs0[n])
    ax.text(1,3.2,r'Scale length of the disk $R_d$',fontsize=12,color='black',horizontalalignment='center',verticalalignment='center')
    ax.text(1,0.6,r'Scale length of the bulge $R_{\rm eff}$',fontsize=12,color='black',horizontalalignment='center',verticalalignment='center')
    ax.set_title('Scale length over time', fontsize = 16)
    ax.set_xlabel('time [Gyr]', fontsize = 14)
    if n==1:
        ax.set_yticklabels([])
    if n==0:
        ax.set_ylabel('Scale length [kpc]', fontsize = 14)
    for k, sim in enumerate(pathlist):
        
        amp1 = np.array([])
        amp2 = np.array([])
        R_eff =[]
        R_d = []
        n = np.array([])
        time = []
 
        simname = glob.glob(sim+'/'+'*.0????')
        #if not os.path.isfile(sim+'_surf_den_test.dat'):
        for name in simname:
            s = pb.load(name)
            s.physical_units()
            new = filt.LowPass('age', s.properties['time'].in_units('Gyr'))
            h1 = s.s[new]
            #h = s.halos()
            #h1 = h[1]
            #h1.physical_units()
            #pb.analysis.angmom.faceon(h1)

            mass_star = (h1.s['mass'].in_units('Msol')).view(np.ndarray)
            rxy_star = np.sqrt((h1.s['x']).view(np.ndarray)**2 + (h1.s['y']).view(np.ndarray)**2)
            srt = np.argsort(rxy_star)
            mass_star = mass_star[srt]
            rxy_star = rxy_star[srt]
        
            rgal = 30 #kpc  #2 * r_90[k]
            eps = h1['eps'].min() #this gives 50 bins
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
            
            popt, pcov = curve_fit(surface_density, R_bins, log_surface_mass_density, bounds=(0, [1e3, 2., 3., 1e3, 5]))
            #print(pcov)
            sersic = models.Sersic1D(popt[0],popt[1],popt[2])#g.amplitude_0, g.r_eff_0, g.n_0) 
            exp1 = models.Exponential1D(popt[3],popt[4])#g.amplitude_1, g.tau_1)
            
            R_eff.append(popt[1])
            R_d.append(popt[4])
            time.append(h1.properties['time'].in_units('Gyr'))


        time2 = np.array(time)
        R_d2 = np.array(R_d)
        R_eff2 = np.array(R_eff)
        srt = np.argsort(time2)
        time3 = time2[srt]
        R_eff3 = R_eff2[srt]
        R_d3 = R_d2[srt]
        

        ax.plot(time3,R_eff3,color=color[k],lw=1, linestyle = '-', label = label[k])
        ax.plot(time3,R_d3,color=color2[k],lw=1)
                #ax.scatter(time3,R_eff3,color=color[k],s=2, linestyle = '-', label = label[k])
                    #ax.scatter(time3,R_d3,color=color2[k],s=2, label = label2[k])
        ax.set_ylim(-0.5,5.5)
        ax.set_xlim(0.1,1.6)
        #ax.set_yticks([0, 1, 2, 3])
        #ax.set_xticks(size = 20)
        #ax.set_yticks(size = 20)
        ax.legend(fontsize=16, loc='center', bbox_to_anchor=(0.5, 0.5))
            #axes[k].scatter(time,R_d,color='black',s=2, label = r'scale length of the disk $R_d$')
                    #axes[k].plot([r_90[k],r_90[k]],[0,4],ls='dashed',color='gray',lw=2)

plt.savefig('scale_length_time_high.pdf', bbox_inches='tight')
