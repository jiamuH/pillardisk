"""Blackbody temperature of the PG 0844+349 non-variable component, and the
thermalizing sphere radius that reproduces it for a 200 Msun embedded star.

The goal here is the SHAPE, not the amplitude. Fit a single blackbody to the
observed constant ("host") component to get its temperature, then ask how big
the thermalized sphere around each embedded star has to be for its surface to
sit at that temperature, given that the star's luminosity is fixed at Eddington:

    sigma T^4 * 4 pi R^2 = L_Edd(M)   ->   R = sqrt(L / (4 pi sigma T^4))

For fixed T the scale factor is linear, so the fit reduces to a one-dimensional
scan over temperature with the amplitude solved in closed form.
"""
import numpy as np
import matplotlib.pyplot as plt

from fluxflux_lib import load_pg0844, SIGMA_SB, L_EDD_PER_MSUN, AU_CM, b_nu, ticks
from pillardisk.pillar_disk import C, ANGSTROM

DL_CM = 287.0 * 3.086e24
NUL = 1e-26 * 4.0 * np.pi * DL_CM**2
MJY = 1e-26                       # erg s^-1 cm^-2 Hz^-1 per mJy
M_STAR = 200.0
N_PILLARS = 300

obs = load_pg0844()
lam, y, e = obs['lam'], obs['host'], obs['host_e']
nu = C / (lam * ANGSTROM)

T_grid = np.linspace(2000.0, 20000.0, 3601)


def fit(weights):
    """Best (T, solid angle) for the given weights. m = B_nu in cgs, y in mJy."""
    chi2, scales = [], []
    for T in T_grid:
        m = b_nu(lam, T)                      # erg/s/cm^2/Hz/ster
        a = np.sum(m * (y * MJY) * weights) / np.sum(m * m * weights)
        chi2.append(np.sum(((a * m - y * MJY) ** 2) * weights))
        scales.append(a)
    i = int(np.argmin(chi2))
    return T_grid[i], scales[i], np.array(chi2)


T_w, om_w, chi2_w = fit(1.0 / (e * MJY) ** 2)          # weighted by the errors
T_u, om_u, chi2_u = fit(1.0 / (y * MJY) ** 2)          # equal fractional weight

print("Blackbody fit to the PG 0844 non-variable component")
print(f"  weighted by quoted errors     : T = {T_w:8.0f} K")
print(f"  equal fractional weighting    : T = {T_u:8.0f} K")
print("  (the quoted errors span four orders of magnitude between UVOT and "
      "LCO, so\n   the weighted fit is dominated by a few LCO points; the "
      "second is more even)")

T_BB = T_u
OMEGA = om_u
print(f"\nAdopting T = {T_BB:.0f} K.")

# emitting area the observed amplitude implies (solid angle x distance^2)
area_proj = OMEGA * DL_CM**2
print(f"Projected emitting area implied by the observed amplitude: "
      f"{area_proj:.3e} cm^2")
print(f"  = a disc of radius {np.sqrt(area_proj/np.pi)/AU_CM:.3e} AU, "
      f"or {np.sqrt(area_proj/np.pi)/1.496e13/2.063e5:.4f} pc")

# sphere radius that puts a 200 Msun star's Eddington luminosity at T_BB
L = M_STAR * L_EDD_PER_MSUN
print(f"\nFor one embedded star of {M_STAR:.0f} Msun, L_Edd = {L:.3e} erg/s.")
print(f"{'geometry':<26}{'area':>14}{'R':>12}")
for label, factor in [("4 pi R^2 (full sphere)", 4.0),
                      ("2 pi R^2 (raised bump)", 2.0),
                      ("pi R^2  (flat patch)", 1.0)]:
    A = L / (SIGMA_SB * T_BB**4)
    R = np.sqrt(A / (factor * np.pi))
    print(f"{label:<26}{A:>14.3e}{R/AU_CM:>11.2f} AU")

# how much of the observed amplitude that provides
R_sph = np.sqrt(L / (4.0 * np.pi * SIGMA_SB * T_BB**4))
f_one = np.pi * b_nu(lam, T_BB) * (R_sph / DL_CM) ** 2 / MJY      # mJy per star
print(f"\nAmplitude check ({N_PILLARS} such stars, shape fixed by T):")
print(f"{'lam':>8}{'model':>11}{'observed':>11}{'ratio':>9}")
for k in range(len(lam)):
    print(f"{lam[k]:>8.0f}{N_PILLARS*f_one[k]:>11.3e}{y[k]:>11.3f}"
          f"{N_PILLARS*f_one[k]/y[k]:>9.4f}")
n_need = np.median(y / f_one)
print(f"  median number of such spheres needed for the amplitude: {n_need:.3e}")

# --- figure: shape comparison, amplitude normalized out --------------------
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

WL = np.logspace(np.log10(1200), np.log10(25000), 300)
nuWL = C / (WL * ANGSTROM)
fig, ax = plt.subplots(figsize=(8.5, 7))
ax.errorbar(lam, nu*y*NUL, yerr=nu*e*NUL, fmt='s', color='darkred', mec='k',
            mew=0.8, ms=9, ecolor='darkred', ls='none', zorder=10,
            label=r'$\rm PG\,0844:~non{-}variable~component$')
for T, om, ls, lab in [(T_u, om_u, '-', rf'$T = {T_u:.0f}~{{\rm K}}$'),
                       (T_w, om_w, '--', rf'$T = {T_w:.0f}~{{\rm K}}$')]:
    mod = om * b_nu(WL, T) / MJY
    ax.plot(WL, nuWL*mod*NUL, ls, color='royalblue', lw=2.5, alpha=0.9,
            label=lab)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_ylim(0.15*(nu*y*NUL).min(), 4.0*(nu*y*NUL).max())
ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=20)
ax.set_ylabel(r'$\nu L_\nu~\rm [erg~s^{-1}]$', fontsize=20)
ax.legend(fontsize=13, frameon=False, loc='lower right')
ticks(ax)
plt.tight_layout()
plt.savefig('plots/fit_pg0844_blackbody.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/fit_pg0844_blackbody.png")
