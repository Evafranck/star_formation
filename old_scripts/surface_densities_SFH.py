import pynbody
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy as np
from pynbody.analysis import profile
import matplotlib.gridspec as gd
from pynbody import units as units
from pynbody import array
import pynbody.filt as f

import warnings


import logging
logger = logging.getLogger('pynbody.plot.stars')

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

s_low_mas = pynbody.load('low_master_iso/low.01000')
s_med_mas = pynbody.load('low_padoan_iso/low.01000')
s_low_sem = pynbody.load('low_semenov_iso/low.01000')
s_med_sem = pynbody.load('low_evans_iso/low.01000')
s_low_mas2 = pynbody.load('med_master_iso/med.01000')
s_med_mas2 = pynbody.load('med_padoan_iso/med.01000')
s_low_sem2 = pynbody.load('med_semenov_iso/med.01000')
s_med_sem2 = pynbody.load('med_evans_iso/med.01000')

s_low_mas.physical_units()
s_med_mas.physical_units()
s_low_sem.physical_units()
s_med_sem.physical_units()
s_low_mas2.physical_units()
s_med_mas2.physical_units()
s_low_sem2.physical_units()
s_med_sem2.physical_units()

lmnew = f.LowPass('age', s_low_mas.properties['time'].in_units('Gyr'))
mmnew = f.LowPass('age', s_med_mas.properties['time'].in_units('Gyr'))
lsnew = f.LowPass('age', s_low_sem.properties['time'].in_units('Gyr'))
msnew = f.LowPass('age', s_med_sem.properties['time'].in_units('Gyr'))
lmnew2 = f.LowPass('age', s_low_mas2.properties['time'].in_units('Gyr'))
mmnew2 = f.LowPass('age', s_med_mas2.properties['time'].in_units('Gyr'))
lsnew2 = f.LowPass('age', s_low_sem2.properties['time'].in_units('Gyr'))
msnew2 = f.LowPass('age', s_med_sem2.properties['time'].in_units('Gyr'))

ps = profile.Profile(s_low_mas.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')
ps2 = profile.Profile(s_med_mas.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')
ps3 = profile.Profile(s_low_sem.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')
ps4 = profile.Profile(s_med_sem.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')
ps_2 = profile.Profile(s_low_mas2.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')
ps22 = profile.Profile(s_med_mas2.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')
ps32 = profile.Profile(s_low_sem2.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')
ps42 = profile.Profile(s_med_sem2.s, nbins=50, rmin='0 kpc', rmax = '50 kpc')


ps_new = pynbody.analysis.profile.Profile(s_low_mas.s[lmnew], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps2_new = pynbody.analysis.profile.Profile(s_med_mas.s[mmnew], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps3_new = pynbody.analysis.profile.Profile(s_low_sem.s[lsnew], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps4_new = pynbody.analysis.profile.Profile(s_med_sem.s[msnew], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps_vert_new = profile.VerticalProfile(s_low_mas.s[lmnew], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)
ps2_vert_new = profile.VerticalProfile(s_med_mas.s[mmnew], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)
ps3_vert_new = profile.VerticalProfile(s_low_sem.s[lsnew], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)
ps4_vert_new = profile.VerticalProfile(s_med_sem.s[msnew], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)

ps_new2 = pynbody.analysis.profile.Profile(s_low_mas2.s[lmnew2], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps2_new2 = pynbody.analysis.profile.Profile(s_med_mas2.s[mmnew2], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps3_new2 = pynbody.analysis.profile.Profile(s_low_sem2.s[lsnew2], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps4_new2 = pynbody.analysis.profile.Profile(s_med_sem2.s[msnew2], nbins=20, rmin='0 kpc', rmax = '20 kpc')
ps_vert_new2 = profile.VerticalProfile(s_low_mas2.s[lmnew2], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)
ps2_vert_new2 = profile.VerticalProfile(s_med_mas2.s[mmnew2], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)
ps3_vert_new2 = profile.VerticalProfile(s_low_sem2.s[lsnew2], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)
ps4_vert_new2 = profile.VerticalProfile(s_med_sem2.s[msnew2], '0.1 kpc', '20 kpc', '5 kpc', ndim = 2, nbins = 20)


rbins = [ps_new['rbins'].in_units('kpc'), ps2_new['rbins'].in_units('kpc'), ps3_new['rbins'].in_units('kpc'),ps4_new['rbins'].in_units('kpc'), ps_new2['rbins'].in_units('kpc'), ps2_new2['rbins'].in_units('kpc'), ps3_new2['rbins'].in_units('kpc'),ps4_new2['rbins'].in_units('kpc')]
zbins = [ps_vert_new['rbins'].in_units('kpc'), ps2_vert_new['rbins'].in_units('kpc'), ps3_vert_new['rbins'].in_units('kpc'),ps4_vert_new['rbins'].in_units('kpc'), ps_vert_new2['rbins'].in_units('kpc'), ps2_vert_new2['rbins'].in_units('kpc'), ps3_vert_new2['rbins'].in_units('kpc'),ps4_vert_new2['rbins'].in_units('kpc')]
surf_den = [np.log10(ps_new['density'].in_units('Msol pc^-2')), np.log10(ps2_new['density'].in_units('Msol pc^-2')), np.log10(ps3_new['density'].in_units('Msol pc^-2')), np.log10(ps4_new['density'].in_units('Msol pc^-2')), np.log10(ps_new2['density'].in_units('Msol pc^-2')), np.log10(ps2_new2['density'].in_units('Msol pc^-2')), np.log10(ps3_new2['density'].in_units('Msol pc^-2')), np.log10(ps4_new2['density'].in_units('Msol pc^-2'))]
vert_surf_den = [np.log10(ps_vert_new['density'].in_units('Msol pc^-2')), np.log10(ps2_vert_new['density'].in_units('Msol pc^-2')), np.log10(ps3_vert_new['density'].in_units('Msol pc^-2')), np.log10(ps4_vert_new['density'].in_units('Msol pc^-2')), np.log10(ps_vert_new2['density'].in_units('Msol pc^-2')), np.log10(ps2_vert_new2['density'].in_units('Msol pc^-2')), np.log10(ps3_vert_new2['density'].in_units('Msol pc^-2')), np.log10(ps4_vert_new2['density'].in_units('Msol pc^-2'))]
colorlist = ['blue','orange', 'green', 'red', 'blue','orange', 'green', 'red']
linestyle = ['-','-', '-','-','--', '--', '--', '--']
title_list = ['SFH', 'Radial surface density profile', 'Vertical surface density profile']
labellist = ['Threshold-based model, Low', 'Padoan (2012), Low', 'Semenov (2016), Low', 'Evans (2022), Low', 'Threshold-based model, Medium', 'Padoan (2012), Medium', 'Semenov (2016), Medium', 'Evans (2022), Medium']
fig = plt.figure(figsize = (17,5))
gs0 = gd.GridSpec(1, 3, figure=fig, width_ratios = [1, 1, 1], height_ratios = [1])

for n in range(3):
    ax = fig.add_subplot(gs0[n])
    if (n==0):
        sfh(s_low_mas, label = 'Threshold-based model, Low', lw = 1, color = 'blue')
        sfh(s_med_mas, label = 'Padoan et al. (2012), Low', lw = 1, color = 'orange')
        sfh(s_low_sem,label = 'Semenov et al. (2016), Low',  lw = 1,color = 'green')
        sfh(s_med_sem, label = 'Evans et al. (2022), Low', lw = 1, color = 'red')
        sfh(s_low_mas2, label = 'Threshold-based model, Medium', lw = 1, ls = '--')
        sfh(s_med_mas2, label = 'Padoan et al. (2012), Medium', lw = 1, ls = '--')
        sfh(s_low_sem2,label = 'Semenov et al. (2016), Medium', lw = 1, ls = '--')
        sfh(s_med_sem2, label = 'Evans et al. (2022), Medium', lw = 1, ls = '--')
        ax.set_title('a) Star formation history', fontsize = 14)
        #ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlim(0, 1.2)
        ax.set_ylim(0,16)
        ax.legend(fontsize = 12)
    if (n==1):
        for i in range(8):
            plt.plot(rbins[i], surf_den[i], lw=1, label = labellist[i], c = colorlist[i], ls = linestyle[i])
        ax.set_xlabel('R [kpc]', fontsize = 14)
        ax.set_xlim(0, 19)
        ax.set_ylim(-4, 3)
        ax.set_ylabel(r'log($\Sigma_{\star}$) [M$_{\odot}$ kpc$^{-2}$]', fontsize = 14)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_title('b) Radial surface density profile', fontsize = 14)
    if (n==2):
        for i in range(8):
            plt.plot(zbins[i], vert_surf_den[i], lw = 1, label = labellist[i], c = colorlist[i], ls = linestyle[i])
        ax.set_xlabel('z [kpc]', fontsize = 14)
        ax.set_xlim(0, 4)
        ax.set_ylim(-3, 3)
        ax.set_ylabel(r'log($\Sigma_{\star}$) [M$_{\odot}$ kpc$^{-2}$]', fontsize = 14)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_title('c) Vertical surface density profile', fontsize = 14)
plt.savefig('SFH_surf_den_iso_all_med_res.pdf', bbox_inches = 'tight')

