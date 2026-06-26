"""Effective blackbody temperature of an AGN accretion-disk hotspot ("pillar")
heated by an embedded star accreting at its Eddington luminosity.

The star's full Eddington luminosity is thermalized and re-emitted from a heated
region of radius R_patch. The emitting area depends on the assumed geometry:
    pi R^2    : flat one-sided patch
    2 pi R^2  : raised 3-D bump radiating into the upper hemisphere  (PRIMARY)
    4 pi R^2  : full sphere (stellar convention)
T = (L / (sigma * A))^(1/4), so at fixed L the temperature scales as R^(-1/2).

Edit the parameters block to try different stellar masses or patch radii.
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.constants import sigma_sb, L_sun, R_sun

# ----------------------------- parameters --------------------------------
M_star  = 200 * u.M_sun      # embedded star mass (paper fiducial)
R_patch = 10 * u.au          # radius of the heated patch
T_sun   = 5772 * u.K         # IAU nominal solar effective temperature
# -------------------------------------------------------------------------

# Eddington luminosity: L_Edd = 1.26e38 (M / M_sun) erg/s
M_ratio = (M_star / u.M_sun).to_value(u.dimensionless_unscaled)
L_edd = 1.26e38 * M_ratio * u.erg / u.s


def temperature(L, R, area_factor):
    """Effective temperature for emitting area A = area_factor * pi R^2."""
    A = area_factor * np.pi * R**2
    flux = (L / A).to(u.erg / u.s / u.cm**2)
    T = ((flux / sigma_sb) ** 0.25).to(u.K)
    return A.to(u.cm**2), flux, T


geometries = [("pi R^2  (flat patch)", 1.0),
              ("2 pi R^2 (raised bump)", 2.0),
              ("4 pi R^2 (full sphere)", 4.0)]

print("INPUTS")
print(f"  M_star  = {M_star:.0f}")
print(f"  R_patch = {R_patch:.1f}  ({R_patch.to(u.cm):.3e})")
print(f"  L_Edd   = {L_edd:.3e}  ({L_edd.to(u.L_sun):.3e})")
print()
print(f"{'geometry':<24}{'area [cm^2]':>14}{'flux [cgs]':>14}{'T [K]':>10}{'T/Tsun':>9}")
for label, factor in geometries:
    A, flux, T = temperature(L_edd, R_patch, factor)
    star = "  <-- primary" if factor == 2.0 else ""
    print(f"{label:<24}{A.value:>14.3e}{flux.value:>14.3e}{T.value:>10.0f}"
          f"{(T/T_sun).to_value(u.dimensionless_unscaled):>9.3f}{star}")

# cross-check the primary (2 pi R^2) against solar scaling:
#   T/T_sun = [ (L/A) / (L_sun/(4 pi R_sun^2)) ]^(1/4)
_, flux2, T2 = temperature(L_edd, R_patch, 2.0)
F_solar = (L_sun / (4 * np.pi * R_sun**2)).to(u.erg / u.s / u.cm**2)
T2_solar = T_sun * ((flux2 / F_solar) ** 0.25)
print(f"\nsolar-scaling check (2 pi R^2): {T2_solar:.0f}  "
      f"(ratio to Stefan-Boltzmann = "
      f"{(T2_solar/T2).to_value(u.dimensionless_unscaled):.5f})")

# size dependence: T ~ R^(-1/2)  (raised-bump, 2 pi R^2)
print("\nSIZE DEPENDENCE (2 pi R^2):  T proportional to R^(-1/2)")
print(f"{'R_patch [AU]':>12}{'T [K]':>10}{'T/Tsun':>9}")
for R_au in [1, 2, 5, 10, 20, 50, 100]:
    _, _, T = temperature(L_edd, R_au * u.au, 2.0)
    print(f"{R_au:>12d}{T.value:>10.0f}"
          f"{(T/T_sun).to_value(u.dimensionless_unscaled):>9.3f}")

# inverse: stellar mass required for a target temperature (hemisphere, R_patch)
# L_req = A sigma T^4 ; M_star = L_req / [1.26e38 erg/s]
print("\nMASS REQUIRED FOR A TARGET T (2 pi R^2, R = 10 AU):")
for T_target in [7000, 8000, 9000] * u.K:
    A = 2 * np.pi * R_patch**2
    L_req = (A * sigma_sb * T_target**4).to(u.erg / u.s)
    M_req = (L_req / (1.26e38 * u.erg / u.s)).to_value(u.dimensionless_unscaled)
    print(f"  T = {T_target:.0f}  ->  M_star = {M_req:.0f} M_sun")

# --- plot T vs R for all three geometries ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

R_grid = np.logspace(np.log10(0.5), np.log10(200), 200) * u.au
colors = {1.0: 'royalblue', 2.0: 'orangered', 4.0: 'seagreen'}
labels = {1.0: r'$\pi R^2~\rm (flat~patch)$',
          2.0: r'$2\pi R^2~\rm (raised~bump)$',
          4.0: r'$4\pi R^2~\rm (sphere)$'}
fig, ax = plt.subplots(figsize=(9, 7))
for factor in (1.0, 2.0, 4.0):
    _, _, T = temperature(L_edd, R_grid, factor)
    lw = 3.5 if factor == 2.0 else 2.5
    ax.plot(R_grid.to_value(u.au), T.to_value(u.K), '-', color=colors[factor],
            lw=lw, label=labels[factor])
ax.axvline(10, color='gray', ls=':', lw=1.5, alpha=0.7)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$\rm pillar~radius~R~[AU]$', fontsize=20)
ax.set_ylabel(r'$\rm hotspot~temperature~T~[K]$', fontsize=20)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=15, frameon=False, title=r'$M_\star=500~M_\odot,~L=L_{\rm Edd}$')
plt.tight_layout()
plt.savefig('plots/hotspot_temperature_vs_R.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/hotspot_temperature_vs_R.png")
