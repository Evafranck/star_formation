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

# Create a list of simulation paths
simulations = ['master', 'semenov', 'evans', 'federrath']

# Create a list of simulation labels and colors
sim_labels = [r'a) Threshold-based model', r'b) Semenov et al. (2016)', r'c) Evans et al. (2022)', r'd) Federrath et al. (2014)']
colorlist = ['blue','orange', 'green', 'red']

# Calculate the Kennicutt-Schmidt law for a given simulation (modified from pynbody documentation)
def schmidtlaw(sim, filename=None, pretime='50 Myr',
			   diskheight='3 kpc', rmax='20 kpc', compare=True,
			   radial=True, bins=10, **kwargs):

	if not radial:
		raise NotImplementedError("Sorry, only radial Schmidt law currently supported")


	if isinstance(pretime, str):
		pretime = units.Unit(pretime)

	# select stuff
	diskgas = sim.gas[f.Disc(rmax, diskheight)]
	diskstars = sim.star[f.Disc(rmax, diskheight)]

	youngstars = np.where(diskstars['tform'].in_units("Myr") >
						  sim.properties['time'].in_units(
							  "Myr", **sim.conversion_context())
						  - pretime.in_units('Myr'))[0]

	# calculate surface densities
	if radial:
		ps = profile.Profile(diskstars[youngstars], nbins=bins)
		pg = profile.Profile(diskgas, nbins=bins)
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
		ysigma = 2.5e-4 * xsigma ** 1.5  # Kennicutt (1998)
		xbigiel = np.logspace(1, 2, 10)
		ybigiel = 10. ** (-2.1) * xbigiel ** 1.0   # Bigiel et al (2007)
	
	return gas_dens, star_dens, xsigma, ysigma, xbigiel, ybigiel
  

# calculate the Kennicutt-Schmidt law for all simulations
def KS(mod):
    s_all = pynbody.load('../high'+'_'+mod+'_iso/' + 'high.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    #disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    #s = s_all[disk]
    s = s_all
    KS_gas_dens.append(schmidtlaw(s, pretime = '100 Myr', compare = True, bins = 20)[0])
    KS_star_dens.append(schmidtlaw(s, pretime = '100 Myr', compare = True, bins = 20)[1])
    KS_xsigma.append(schmidtlaw(s, pretime = '100 Myr', compare = True, bins = 20)[2])
    KS_ysigma.append(schmidtlaw(s, pretime = '100 Myr', compare = True, bins = 20)[3])

for sim_path in simulations:
    KS(sim_path)

# plot the Kennicutt-Schmidt law for all simulations
fig = plt.figure(figsize=(10, 10))
grid = gd.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

for n in range(4):
    ax = plt.subplot(grid[n])
    ax.loglog(KS_gas_dens[n], KS_star_dens[n], "+", label = 'simulation data', color = colorlist[n])
    ax.loglog(KS_xsigma[n], KS_ysigma[n], label='Kennicutt (1998)')
    ax.vlines(1, 1.5, 7.9, ls = 'dashed', color = 'grey', linewidths = 1, label = r'0.1% $\epsilon_{\rm{ff}}$')
    ax.set_xlabel('$\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]')
    ax.set_ylabel('$\Sigma_{SFR}$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]')
    ax.legend(loc = 'lower right')
    ax.set_title(sim_labels[n], fontsize = 16)
plt.tight_layout()
plt.show()
plt.savefig('KS_law_separate.pdf', bbox_inches='tight')

#ax.loglog(KS_xbigiel[n], KS_ybigiel[n], linestyle="dashed", label='Bigiel et al (2007)')
#ax.vlines(1, 1.5, 7.9, ls = 'dashed', color = 'grey', linewidths = 1, label = r'1% $\epsilon_{\rm{ff}}$')
#ax.vlines(1, 1.5, 7.9, ls = 'dashed', color = 'grey', linewidths = 1, label = r'10% $\epsilon_{\rm{ff}}$')