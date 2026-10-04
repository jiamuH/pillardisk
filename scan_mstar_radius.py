"""Which (M_star, R) combinations put the thermalized sphere at roughly the
temperature of the PG 0844+349 non-variable component?

The embedded star radiates at Eddington, so L is proportional to M, and the
sphere surface sits at sigma T^4 = L / A. For A = f pi R^2 that gives

    T = [ 1.26e38 (M/Msun) / (f pi sigma R^2) ]^(1/4)  ->  T ~ M^(1/4) R^(-1/2)

so only the combination M / R^2 is constrained. The observed temperature cannot
pin the stellar mass on its own; it fixes a one-parameter family, and this maps
it. The target band is taken from the two weightings of the blackbody fit to the
observed component, which bracket the plausible temperature.
"""
import numpy as np
import matplotlib.pyplot as plt

from fluxflux_lib import load_pg0844, SIGMA_SB, L_EDD_PER_MSUN, AU_CM, b_nu, ticks
from pillardisk.pillar_disk import C, ANGSTROM

MJY = 1e-26
GEOM = {'4 pi R^2 (sphere)': 4.0,
        '2 pi R^2 (raised bump)': 2.0,
        'pi R^2 (flat patch)': 1.0}
GEOM_PRIMARY = 2.0                 # the raised-bump convention used in the paper

# --- temperature band from the blackbody fit -------------------------------
obs = load_pg0844()
lam, y, e = obs['lam'], obs['host'], obs['host_e']
T_grid = np.linspace(2000.0, 20000.0, 3601)


def fit_T(weights):
    chi2 = []
    for T in T_grid:
        m = b_nu(lam, T)
        a = np.sum(m * (y * MJY) * weights) / np.sum(m * m * weights)
        chi2.append(np.sum(((a * m - y * MJY) ** 2) * weights))
    return T_grid[int(np.argmin(chi2))]


T_LO = fit_T(1.0 / (e * MJY) ** 2)          # error weighted
T_HI = fit_T(1.0 / (y * MJY) ** 2)          # equal fractional weight
print(f"Blackbody fit to the observed non-variable component: "
      f"T = {T_LO:.0f}-{T_HI:.0f} K\n"
      f"  (the two ends are the error-weighted and equal-fractional fits; the "
      f"spread\n   is the honest uncertainty, so any (M, R) landing in this "
      f"band is acceptable)\n")


def radius_for(m_star, T, factor=GEOM_PRIMARY):
    """Sphere radius (cm) giving surface temperature T for an Eddington star."""
    return np.sqrt(m_star * L_EDD_PER_MSUN / (factor * np.pi * SIGMA_SB * T**4))


def temperature_for(m_star, r_au, factor=GEOM_PRIMARY):
    A = factor * np.pi * (r_au * AU_CM) ** 2
    return (m_star * L_EDD_PER_MSUN / (A * SIGMA_SB)) ** 0.25


print("Scaling relation (raised-bump geometry):")
T_ref = temperature_for(200.0, 10.0)
print(f"  T = {T_ref:.0f} K x (M/200 Msun)^(1/4) x (R/10 AU)^(-1/2)\n")

print("Radius that lands in the band, per stellar mass (raised bump, 2 pi R^2):")
print(f"{'M_star':>9}{'R at ' + f'{T_HI:.0f} K':>14}{'R at ' + f'{T_LO:.0f} K':>14}"
      f"{'R range':>16}")
for m in [50.0, 100.0, 200.0, 300.0, 500.0, 1000.0]:
    r_hi = radius_for(m, T_HI) / AU_CM
    r_lo = radius_for(m, T_LO) / AU_CM
    print(f"{m:>9.0f}{r_hi:>13.2f} AU{r_lo:>13.2f} AU"
          f"{f'{r_hi:.1f}-{r_lo:.1f}':>16}")

print("\nSame, for the other emitting-area conventions at M = 200 Msun:")
for label, f in GEOM.items():
    print(f"  {label:<26}"
          f"{radius_for(200.0, T_HI, f)/AU_CM:>7.2f}-"
          f"{radius_for(200.0, T_LO, f)/AU_CM:<7.2f} AU")

print("\nOnly M/R^2 is constrained, so the mass is not determined by the "
      "temperature:\nany point on the band is equally good, and picking M "
      "requires a separate argument.")

# --- figure -----------------------------------------------------------------
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

M = np.logspace(np.log10(20), np.log10(2000), 400)
fig, ax = plt.subplots(figsize=(8.5, 7))
ax.fill_between(M, radius_for(M, T_HI)/AU_CM, radius_for(M, T_LO)/AU_CM,
                color='royalblue', alpha=0.30, lw=0,
                label=rf'$2\pi R^2~{{\rm (raised~bump)}},~'
                      rf'T = {T_LO:.0f}-{T_HI:.0f}~{{\rm K}}$')
ax.plot(M, radius_for(M, 0.5*(T_LO+T_HI))/AU_CM, '-', color='royalblue',
        lw=2.5)
for (f, tex), ls in [((4.0, r'$4\pi R^2~\rm (full~sphere)$'), '--'),
                     ((1.0, r'$\pi R^2~\rm (flat~patch)$'), ':')]:
    ax.plot(M, radius_for(M, 0.5*(T_LO+T_HI), f)/AU_CM, ls, color='dimgray',
            lw=2.0, label=tex)
ax.plot([200.0], [10.0], marker='*', color='orangered', ms=22, mec='k',
        mew=1.0, ls='none', zorder=10,
        label=r'$\rm paper~fiducial:~200~M_\odot,~10~AU$')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$M_\star~\rm [M_\odot]$', fontsize=20)
ax.set_ylabel(r'$\rm thermalized~sphere~radius~R~[AU]$', fontsize=20)
ax.legend(fontsize=13, frameon=False, loc='upper left')
ticks(ax)
plt.tight_layout()
plt.savefig('plots/scan_mstar_radius.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/scan_mstar_radius.png")
