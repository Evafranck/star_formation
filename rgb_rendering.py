import pynbody
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gd
import pynbody.filt as f


rgb_arrays_faceon = []
rgb_arrays_sideon = []

# Load a slice of the simulation snapshots faceon and sideon
def load_sim_faceon(mod):
    s_all = pynbody.load('../' + mod + '/halo.00128')
    pynbody.analysis.angmom.faceon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s = s_all[disk]
    rgb_arrays_faceon.append(pynbody.plot.stars.render(s_all, width='30 kpc', resolution=500, ret_im=True))

def load_sim_sideon(mod):
    s_all = pynbody.load('../' +mod + '/halo.00128')
    pynbody.analysis.angmom.sideon(s_all)
    s_all.physical_units()
    disk = f.LowPass('r', '30 kpc') & f.BandPass('z', '-5 kpc', '5 kpc')
    s = s_all[disk]
    rgb_arrays_sideon.append(pynbody.plot.stars.render(s_all, width='30 kpc', resolution=500, ret_im=True))

#simulations = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan']

#sim_labels = ['Threshold-based model', 'Federrath & Klessen (2012)', 'Hopkins et al. (2013) with' + '\n' + 'efficiency of Padoan et al. (2012)', 'Hopkins et al. (2013) with' + '\n' + r'$\alpha_{\mathrm{vir}}$ threshold', r'Hopkins et al. (2013) with' + '\n' + r' $\alpha_{\mathrm{vir}}$ of Padoan et al. (2012)']

#simulations = ['threshold', 'federrath', 'hopkins', 'hopkins_alpha', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
simulations = ['threshold_alpha008','semenov_1e6_alpha008', 'semenov_alpha008', 'semenov_cstar_cut', 'federrath_1e6_alpha008']
#, 'federrath_alpha008', 'federrath_cstar_cut'] 
#simulations = ['threshold_alpha008', 'threshold_1e6_alpha008', 'hopkins_alpha_padoan_alpha008', 'hopkins_alpha008', 'hopkins_alpha_padoan', 'hopkins_alpha_alpha008']
sim_labels = simulations
for sim_path in simulations:
    load_sim_faceon(sim_path)
for sim_path in simulations:    
    load_sim_sideon(sim_path)

fig = plt.figure(figsize = (12.5, 4))
gs0 = gd.GridSpec(2, 5)
#gs0.update(hspace=0.00, wspace=0.00)

for n in range(5):
    ax1 = plt.subplot(gs0[n])
    
    # Display the RGB rendering in the subplot
    ax1.imshow(rgb_arrays_faceon[n], origin='lower')
    ax1.axis('off')
    
    ax1.set_title(sim_labels[n], fontsize = 12)

for n in range(5, 10):
    ax1 = plt.subplot(gs0[n])
    
    # Display the RGB rendering in the subplot
    ax1.imshow(rgb_arrays_sideon[n-5], origin='lower')
    ax1.axis('off')
    ax1.set_ylim(0,200)
plt.tight_layout()
plt.show()
plt.savefig('rgb_rendering.pdf', bbox_inches='tight')
