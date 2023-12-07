import pynbody
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gd
import pynbody.filt as f
import glob
from pynbody.analysis import profile
from pynbody import units as units
from pynbody import array
import pynbody.plot as pp

KS_gas_dens = []
KS_star_dens = []
KS_xsigma = []
KS_ysigma = []
KS_xbigiel = []
KS_ybigiel = []

# calculate SFR with t_ff and epsilon_ff = 0.1%, 1%, 10% 	
x = np.logspace(-1.5, 1, 100)  # Generates 100 points from 10^0 to 10^2 Msol/pc^2
G = 6.7*10**(-11)*units.m**3/(units.kg*units.s**2)
t_ff = 1e8 #np.sqrt(3*np.pi/(32*G.in_units('pc**3 Msol**-1 yr**-2')*x)) # in yr
SFR = x*10**6/t_ff # in Msol/yr/kpc^2


# Create a list of simulation paths
simulations = ['master', 'semenov', 'evans', 'federrath']

# Create a list of simulation labels and colors
sim_labels = [r'Threshold-based model', r'Semenov et al. (2016)', r'Evans et al. (2022)', r'Federrath et al. (2012)']
colorlist = ['blue','orange', 'green', 'red']

# Calculate the Kennicutt-Schmidt law for a given simulation (modified from pynbody documentation)
def schmidtlaw(sim, filename=None, pretime='10 Myr',
			   diskheight='2 kpc', rmax='30 kpc', compare=True,
			   radial=True, bins=10, **kwargs):

	if not radial:
		raise NotImplementedError("Sorry, only radial Schmidt law currently supported")


	if isinstance(pretime, str):
		pretime = units.Unit(pretime)
  
	gas = sim.gas[f.Disc(rmax, diskheight)]
	cold = f.LowPass('temp', '30000 K')
	diskgas = gas[cold]
	diskstars = sim.star[f.Disc(rmax, diskheight)]
    
	youngstars = np.where(diskstars['tform'].in_units("Myr") >
						  sim.properties['time'].in_units(
							  "Myr", **sim.conversion_context())
						  - pretime.in_units('Myr'))[0]

	# calculate surface densities
	if radial:
		ps = profile.Profile(diskstars[youngstars], bins=np.linspace(0, 30, 60))
		pg = profile.Profile(diskgas, bins=np.linspace(0, 30, 60))

	else:
		# make bins 2 kpc
		nbins = rmax * 2 / binsize
		pg, x, y = np.histogram2d(diskgas['x'], diskgas['y'], bins=nbins,
								  weights=diskgas['mass'],
								  range=[(-rmax, rmax), (-rmax, rmax)])
		ps, x, y = np.histogram2d(diskstars[youngstars]['x'],
								  diskstars[youngstars]['y'],
								  weights=diskstars['mass'],
								  bins=nbins, range=[(-rmax, rmax), (-rmax, rmax)])

	gas_dens = pg['density'].in_units('Msol pc^-2')
	star_dens = ps['density'].in_units('Msol kpc^-2')/pretime/1e6

 

	if compare:
		xsigma = np.logspace(np.log10(pg['density'].in_units('Msol pc^-2').min()),
							 np.log10(
								 pg['density'].in_units('Msol pc^-2').max()),
							 100)
		ysigma = 2.5e-4 * xsigma ** 1.4  # Kennicutt (1998)
		xbigiel = np.logspace(0, 1, 10)
		ybigiel = 10. ** (-2.1) * xbigiel ** 1.0   # Bigiel et al (2007)
	
	return gas_dens, star_dens, xsigma, ysigma, xbigiel, ybigiel
  

# calculate the Kennicutt-Schmidt law for all simulations
def KS(mod):
    s_all = pynbody.load('../high'+'_'+mod+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    s = s_all[disk]
    KS_gas_dens.append(schmidtlaw(s, pretime = '10 Myr', compare = True)[0])
    KS_star_dens.append(schmidtlaw(s, pretime = '10 Myr', compare = True)[1])
    KS_xsigma.append(schmidtlaw(s, pretime = '10 Myr', compare = True)[2])
    KS_ysigma.append(schmidtlaw(s, pretime = '10 Myr', compare = True)[3])
    KS_xbigiel.append(schmidtlaw(s, pretime = '10 Myr', compare = True)[4])
    KS_ybigiel.append(schmidtlaw(s, pretime = '10 Myr', compare = True)[5])

for sim_path in simulations:
    KS(sim_path)
    
# plot the Kennicutt-Schmidt law for all simulations
fig = plt.figure(figsize=(10, 10))

for n in range(len(simulations)):
    plt.loglog(KS_gas_dens[n], KS_star_dens[n], "+", label = sim_labels[n], color = colorlist[n])
plt.loglog(KS_xsigma[0], KS_ysigma[0], label='Kennicutt (1998)')
plt.loglog(x, 0.1*SFR, ls = '-', color = 'grey', label = r'10% $\epsilon_{\rm{ff}}$')
plt.loglog(x, 0.01*SFR, ls = '--', color = 'grey', label = r'1% $\epsilon_{\rm{ff}}$')
plt.loglog(x, 0.001*SFR, ls = ':', color = 'grey', label = r'0.1% $\epsilon_{\rm{ff}}$')
plt.xlabel('$\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]')
plt.ylabel('$\Sigma_{SFR}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]')
plt.legend(loc = 'lower right', fontsize = 14)
plt.title('Kennicutt-Schmidt-relation', fontsize = 16)
plt.tight_layout()
plt.show()
plt.savefig('KS_law_high.pdf', bbox_inches='tight')

