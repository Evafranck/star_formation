import pynbody
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
from pynbody.analysis import profile
import matplotlib.gridspec as gd
from pynbody import units as units
from pynbody import array
import pynbody.filt as f

def sfh(sim, filename=None, massform=True, clear=False, legend=False,
		subplot=False, trange=False, bins=200, **kwargs):
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


#labellist = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with ' + '\n' + r'$\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)']
#colorlist = ['blue','orange', 'green', 'red', 'purple']
#model = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']

model_name = '1e6' # options: 'hopkins', 'semenov', 'federrath_padoan', '1e6'
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

    
labellist = model


plt.figure(figsize = (10,10))
plt.title('Star formation history', fontsize = 16)
plt.xlabel('Time [Gyr]', fontsize = 14)
plt.ylabel('log(SFR) [M$_\odot$ yr$^{-1}$]', fontsize = 14)
#plt.xlim(0, 2.5)
plt.ylim(2, 4e1)
plt.yscale('log')
for n in range(len(model)):
    s = pynbody.load('../'+model[n]+'/halo.00128')
    pynbody.analysis.angmom.faceon(s)
    s.physical_units()
    sfh(s, label = labellist[n], lw = 1) #, color = colorlist[n])
plt.legend(fontsize = 14)
plt.savefig('SFH_time_' + model_name + '.pdf', bbox_inches = 'tight')

