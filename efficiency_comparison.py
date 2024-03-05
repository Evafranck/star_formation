import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf


Mach_range = [0.1, 1.0, 10.0, 100.0]
#label_list = ['Federrath & Klessen (2012) with Mach number = ' + str(Mach), 'Padoan et al. (2012)',  'Semenov et al. (2016)', 'Evans et al. (2022)']
color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

def Evans(alpha):
    return 0.3*np.exp(-2.02*alpha**(1/2))

def Semenov(alpha):
    return 0.9*np.exp(-1.6*alpha**(1/2))

def Padoan(alpha):
    return 0.5*np.exp(-1.6*alpha**(1/2))

def Hopkins(alpha):
    return np.where(alpha < 1, 1, 0)

def Federrath(alpha):
    Mach = np.sqrt(alpha)
    print(Mach.max())
    return 0.09/(2*0.49)*np.exp(3/8*np.log(1+0.4**2*Mach**2))*(1+erf((np.log(1+0.4**2*Mach**2)-np.log(np.pi**2/5*0.19**2*alpha*Mach**2))/np.sqrt(2*np.log(1+0.4**2*Mach**2))))

# Federrath in dependence of Mach
def Federrath_mach(alpha, Mach):
    return 0.09/(2*0.49)*np.exp(3/8*np.log(1+0.4**2*Mach**2))*(1+erf((np.log(1+0.4**2*Mach**2)-np.log(np.pi**2/5*0.19**2*alpha*Mach**2))/np.sqrt(2*np.log(1+0.4**2*Mach**2))))

fig = plt.figure(figsize = (10,8))
alpha = np.logspace(-3, 4, 100)
plt.xscale('log')
plt.ylim(-0.05, 1.05)
for Mach in Mach_range:
    plt.plot(alpha, Federrath_mach(alpha, Mach), label='Federrath & Klessen (2012)' + '\n' + 'with Mach number = ' + str(Mach))    
plt.plot(alpha, Federrath(alpha), label='Federrath & Klessen (2012)')
plt.title('Star Formation Efficiency Comparison in different models', fontsize = 14)
plt.plot(alpha, Padoan(alpha), label='Padoan et al. (2012)')
plt.plot(alpha, Hopkins(alpha), label='Hopkins et al. (2013)')
plt.plot(alpha, Semenov(alpha), label='Semenov et al. (2016)')
plt.plot(alpha, Evans(alpha), label='Evans et al. (2022)')
plt.xlabel(r'Virial Parameter $\alpha$', fontsize = 14)
#fig.tight_layout()
plt.ylabel(r'Star Formation Efficiency $\epsilon_{\rm ff}$', fontsize = 14)
plt.legend(fontsize = 11)
plt.savefig('efficiency_comparison.pdf')

