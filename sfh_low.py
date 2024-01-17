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

def sfh(sim, filename=None, massform=True, clear=False, legend=False,
		subplot=False, trange=False, bins=100, **kwargs):
	'''
	star formation history

	**Optional keyword arguments:**

	   *trange*: list, array, or tuple
		 size(t_range) must be 2. Specifies the time range.

	   *bins*: int
		 number of bins to use for the SFH

	   *massform*: bool
		 decides whether to use original star mass (massform) or final star mass

	   *subplot*: subplot object
		 where to plot SFH

	   *legend*: boolean
		 whether to draw a legend or not

	   *clear*: boolean
		 if False (default), plot on the current axes. Otherwise, clear the figure first.

	By default, sfh will use the formation mass of the star.  In tipsy, this will be
	taken from the starlog file.  Set massform=False if you want the final (observed)
	star formation history

	**Usage:**

	>>> import pynbody.plot as pp
	>>> pp.sfh(s,linestyle='dashed',color='k')


	'''
	import matplotlib.pyplot as pyplot

	if subplot:
		plt = subplot
	else:
		plt = pyplot

	if "nbins" in kwargs:
		bins = kwargs['nbins']

	if 'nbins' in kwargs:
		bins = kwargs['nbins']
		del kwargs['nbins']

	if ((len(sim.g)>0) | (len(sim.d)>0)): simstars = sim.star
	else: simstars = sim

	if trange:
		assert len(trange) == 2
	else:
		trange = [simstars['tform'].in_units(
			"Gyr").min(), simstars['tform'].in_units("Gyr").max()]
	binnorm = 1e-9 * bins / (trange[1] - trange[0])

	trangefilt = f.And(f.HighPass('tform', str(trange[0]) + ' Gyr'),
						  f.LowPass('tform', str(trange[1]) + ' Gyr'))
	tforms = simstars[trangefilt]['tform'].in_units('Gyr')

	if massform:
		try:
			weight = simstars[trangefilt][
				'massform'].in_units('Msol') * binnorm
		except (KeyError, units.UnitsException):
			warnings.warn(
				"Could not load massform array -- falling back to current stellar masses", RuntimeWarning)
			weight = simstars[trangefilt]['mass'].in_units('Msol') * binnorm
	else:
		weight = simstars[trangefilt]['mass'].in_units('Msol') * binnorm

	if clear:
		plt.clf()
	sfhist, thebins, patches = plt.hist(tforms, weights=weight, bins=bins,
										histtype='step', **kwargs)
	if not subplot:
		# don't set the limits
		#plt.ylim(0.0, 1.2 * np.max(sfhist))
		plt.xlabel('Time [Gyr]', fontsize=14)
		plt.ylabel('SFR [M$_\odot$ yr$^{-1}$]', fontsize=14)
	else:
		plt.set_ylim(0.0, 1.2 * np.max(sfhist))

	# Make both axes have the same start and end point.
	if subplot:
		x0, x1 = plt.get_xlim()
	else:
		x0, x1 = plt.gca().get_xlim()

	if legend:
		plt.legend(loc=1)
	if filename:
		logger.info("Saving %s", filename)
		plt.savefig(filename)

	return array.SimArray(sfhist, "Msol yr**-1"), array.SimArray(thebins, "Gyr")

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

sfh_list = []

# Create a list of simulation paths
simulations = ['../threshold', '../federrath', 'master', 'federrath']

# Load a slice of the simulation snapshots faceon and sideon
def load_sim_faceon(mod):
    if (mod == '../threshold' or mod == '../federrath'):
        s_all = pynbody.load(mod + '/halo.00128')
    else:
        s_all = pynbody.load('../low'+'_'+mod+'_iso/' +'low.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    #new = f.LowPass('age', s_all.properties['time'].in_units('Gyr'))
    new = f.LowPass('age', '1.49 Gyr')
    s_new = s_all.s[new]
    cold = f.LowPass('temp', '30000 K') # nur kaltes gas
    s_gas = s_all.g[cold]
    sfh_list.append(s_new)

    
for m in simulations:
    load_sim_faceon(m)
    

#titlelist = ['a) Threshold-based model', 'b) Padoan et al. (2012)', 'c) Federrath & Klessen (2012)' + '\n' + 'with temperature cut', 'd) Federrath & Klessen (2012)' + '\n' + 'without temperature cut', '', '', '', '']
titlelist = ['Threshold-based model' + '\n' + 'Jakob Herpichs ICs', 'Federrath & Klessen (2012)' + '\n' + 'Jakob Herpichs ICs', 'Threshold-based model' + '\n' + 'AGORA ICs', 'Federrath & Klessen (2012)' + '\n' + 'AGORA ICs', '', '', '', ''] # r'Hopkins et al. (2013)' + '\n' + 'with temperature cut', r'Hopkins et al. (2013)' + '\n' + 'without temperature cut'] # 'Federrath & Klessen (2012)' + '\n' + 'with temperature cut']

fig = plt.figure(figsize = (11,5))
gs0 = gd.GridSpec(1, 2, figure=fig, width_ratios = [1, 1])

for n in range(2):
    ax = fig.add_subplot(gs0[n])
    if (n==0):
        for i in range(2):
            sfh(sfh_list[i], lw = 1, label = titlelist[i])
        ax.set_title('a) Star formation history with Jakob Herpichs ICs', fontsize = 12)
        #ax.set_aspect(1./ax.get_data_ratio())
        #ax.set_xlim(0, 1.2)
        #ax.set_ylim(0,16)
        ax.legend(fontsize = 11, loc = 'lower center')
        ax.set_ylabel(r'SFR [M$_{\odot} yr^{-1}$]', fontsize = 14)
        ax.set_xlabel('Time [Gyr]', fontsize = 14)
        
    if (n==1):
        for i in range(2,4):
            sfh(sfh_list[i], lw = 1, label = titlelist[i])
        ax.set_title('b) Star formation history with AGORA ICs', fontsize = 12)
        #ax.set_aspect(1./ax.get_data_ratio())
        #ax.set_xlim(0, 1.2)
        #ax.set_ylim(0,16)
        ax.legend(fontsize = 11, loc = 'upper right')
        ax.set_ylabel(r'SFR [M$_{\odot} yr^{-1}$]', fontsize = 14)
        ax.set_xlabel('Time [Gyr]', fontsize = 14)
        
plt.savefig('sfh_low_ICs_comparison.pdf', bbox_inches = 'tight')



