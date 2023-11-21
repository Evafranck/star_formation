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

rgb_arrays_faceon = []
rgb_arrays_sideon = []

# Load a slice of the simulation snapshots faceon and sideon
def load_sim_faceon(mod):
    s_all = pynbody.load('../high'+'_'+mod+'_iso/' +'high.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    s = s_all[disk]
    rgb_arrays_faceon.append(pynbody.plot.stars.render(s_all, width='30 kpc', resolution=500, ret_im=True))

def load_sim_sideon(mod):
    s_all = pynbody.load('.high'+'_'+mod+'_iso/' + 'low.01000')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    s = s_all[disk]
    rgb_arrays_sideon.append(pynbody.plot.stars.render(s_all, width='30 kpc', resolution=500, ret_im=True))

# Create a list of simulation paths
simulations = ['master', 'semenov', 'evans', 'federrath']

# Create a list of simulation labels (for titles)
sim_labels = [r'a) Threshold-based model', r'b) Semenov et al. (2016)', r'c) Evans et al. (2022)', r'd) Federrath et al. (2014)']

for sim_path in simulations:
    load_sim_faceon(sim_path)
for sim_path in simulations:    
    load_sim_sideon(sim_path)

fig = plt.figure(figsize = (12, 4))
gs0 = gd.GridSpec(2, 4, height_ratios = [1, 1], width_ratios = [1, 1, 1, 1])
gs0.update(hspace=0.00, wspace=0.00)

for n in range(4):
    ax1 = plt.subplot(gs0[0, n])
    ax2 = plt.subplot(gs0[1, n])
    
    # Display the RGB rendering in the subplot
    ax1.imshow(rgb_arrays_faceon[n], origin='lower')
    ax1.axis('off')
    
    ax2.imshow(rgb_arrays_sideon[n], origin='lower')
    ax2.axis('off')
    
    ax1.set_title(sim_labels[n])

plt.tight_layout()
plt.show()
plt.savefig('rgb_rendering_high.pdf', bbox_inches='tight')
