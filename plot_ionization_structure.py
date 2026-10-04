"""Referee reply, Section 4 bullet 7: ionization structure of the illuminated
grid point (log Phi_H = 19, log n_H = 10, solar metallicity, Nagao SED).

Bottom axis: depth into the cloud in cm. Top axis: hydrogen column density
N_H = n_H * depth, with n_H = 10^10 cm^-3 constant through the model.

Needs the Cloudy run ionstruct_n10_phi19 (deck in the same directory).

Run:  python3 plot_ionization_structure.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt

D = '/Users/jiamuh/c23.01/my_models/sed_test_ngc5548'
STEM = 'ionstruct_n10_phi19'
DENSITY = 1e10  # cm^-3
PLOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def read_element(name):
    path = os.path.join(D, '%s_%s.ele' % (STEM, name))
    with open(path) as f:
        header = f.readline().lstrip('#').rstrip('\n').split('\t')
    data = np.loadtxt(path, comments='#')
    return {h.strip(): data[:, i] for i, h in enumerate(header)}


def read_overview():
    path = os.path.join(D, '%s_blr.ovr' % STEM)
    with open(path) as f:
        header = f.readline().lstrip('#').rstrip('\n').split('\t')
    data = np.loadtxt(path, comments='#')
    return {h.strip(): data[:, i] for i, h in enumerate(header)}


hydrogen = read_element('hydrogen')
carbon = read_element('carbon')
magnesium = read_element('magnesium')
overview = read_overview()
depth = hydrogen['depth']

curves = [(hydrogen['H+'], r'$\rm H^{+}$', 'black', '-'),
          (hydrogen['H'], r'$\rm H^{0}$', 'gray', '--'),
          (overview['HeII'], r'$\rm He^{+}$', 'forestgreen', '-'),
          (carbon['C+3'], r'$\rm C^{3+}$', 'orangered', '-'),
          (magnesium['Mg+'], r'$\rm Mg^{+}$', 'royalblue', '-')]

from matplotlib.ticker import MultipleLocator

log_depth = np.log10(depth)

fig, ax = plt.subplots(figsize=(9, 6.5))
for values, label, color, style in curves:
    ax.plot(log_depth, values, color=color, ls=style, lw=3, alpha=0.8, label=label)

ax.set_yscale('log')
ax.set_xlim(log_depth[0], log_depth[-1])
ax.set_ylim(1e-5, 2.)
ax.set_xlabel(r'$\rm log~Depth~[cm]$')
ax.set_ylabel(r'$\rm Ion~fraction$')
ax.xaxis.set_major_locator(MultipleLocator(2.))
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
ax.tick_params(which='both', direction='in', top=False, right=True)
ax.tick_params(which='major', length=8, width=1.5)
ax.tick_params(which='minor', length=4, width=1)
ax.legend(fontsize=13, loc='center left')

top = ax.twiny()
top.set_xlim(log_depth[0] + np.log10(DENSITY), log_depth[-1] + np.log10(DENSITY))
top.set_xlabel(r'$\rm log$ $N_{\rm H}~[\rm cm^{-2}]$')
top.xaxis.set_major_locator(MultipleLocator(2.))
top.xaxis.set_minor_locator(MultipleLocator(0.5))
top.tick_params(which='both', direction='in')
top.tick_params(which='major', length=8, width=1.5)
top.tick_params(which='minor', length=4, width=1)

fig.tight_layout()
out = os.path.join(PLOT_DIR, 'ionization_structure_phi19.png')
fig.savefig(out, dpi=200)
print('saved', out)

for target, label in ((0.5, 'H+ drops below 0.5'), (0.1, 'H+ drops below 0.1')):
    k = np.argmax(hydrogen['H+'] < target)
    if k:
        print('%-22s at depth %.2e cm, log N_H = %.2f'
              % (label, depth[k], np.log10(depth[k] * DENSITY)))
for values, label in ((carbon['C+3'], 'C+3'), (magnesium['Mg+'], 'Mg+')):
    k = int(np.argmax(values))
    print('%-4s peaks at depth %.2e cm, log N_H = %.2f, fraction %.2f'
          % (label, depth[k], np.log10(depth[k] * DENSITY), values[k]))
