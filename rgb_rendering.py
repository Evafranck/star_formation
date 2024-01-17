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
    if (mod == '../threshold' or mod == '../federrath'):
        s_all = pynbody.load(mod + '/halo.00128')
    else:
        s_all = pynbody.load('../med'+'_'+mod+'_iso/' +'med.01000')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    s = s_all[disk]
    rgb_arrays_faceon.append(pynbody.plot.stars.render(s_all, width='20 kpc', resolution=500, ret_im=True))

def load_sim_sideon(mod):
    if (mod == '../threshold' or mod == '../federrath'):
        s_all = pynbody.load(mod + '/halo.00128')
    else:
        s_all = pynbody.load('../med'+'_'+mod+'_iso/' +'med.01000')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-10 kpc', '10 kpc')
    s = s_all[disk]
    #rgb_arrays_sideon.append(pynbody.plot.stars.render(s_all, width='30 kpc', resolution=500, ret_im=True))

# Create a list of simulation paths
simulations = ['../threshold', '../federrath', 'master', 'federrath']

# Create a list of simulation labels (for titles)
sim_labels = ['Threshold-based model' + '\n' + 'Jakob Herpichs ICs', 'Federrath & Klessen (2012)' + '\n' + 'Jakob Herpichs ICs', 'Threshold-based model' + '\n' + 'AGORA ICs', 'Federrath & Klessen (2012)' + '\n' + 'AGORA ICs'] # r'Hopkins et al. (2013)' + '\n' + 'with temperature cut', r'Hopkins et al. (2013)' + '\n' + 'without temperature cut'] # 'Federrath & Klessen (2012)' + '\n' + 'with temperature cut']

for sim_path in simulations:
    load_sim_faceon(sim_path)
for sim_path in simulations:    
    load_sim_sideon(sim_path)

fig = plt.figure(figsize = (12, 12))
gs0 = gd.GridSpec(2, 2, height_ratios = [1, 1], width_ratios = [1, 1])
#gs0.update(hspace=0.00, wspace=0.00)

for n in range(4):
    ax1 = plt.subplot(gs0[n])
    
    # Display the RGB rendering in the subplot
    ax1.imshow(rgb_arrays_faceon[n], origin='lower')
    ax1.axis('off')
    
    ax1.set_title(sim_labels[n], fontsize = 16)

plt.tight_layout()
plt.show()
plt.savefig('rgb_rendering.pdf', bbox_inches='tight')
