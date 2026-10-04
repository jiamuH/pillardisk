"""Can shrinking the thermalizing sphere match the PG 0844 non-variable component?

Shrinking the emitting area raises the temperature at FIXED luminosity, since
sigma T^4 A = L. The question is whether the resulting hotter, smaller emitter
can put more flux into the observed bands. For a fixed total luminosity there is
a hard ceiling: at any wavelength, nu L_nu of a blackbody is at most a fraction
of L_bol, maximised when the blackbody peaks at that wavelength. This scans
temperature (equivalently sphere radius) at fixed L and compares the best
achievable nu L_nu against the observed constant component.
"""
import numpy as np

from fluxflux_lib import load_pg0844, SIGMA_SB, L_EDD_PER_MSUN, AU_CM, b_nu
from pillardisk.pillar_disk import C, ANGSTROM

DL_CM = 287.0 * 3.086e24
NUL = 1e-26 * 4.0 * np.pi * DL_CM**2

N_PILLARS = 300
M_STAR = 200.0
LAM_V = 5468.0

L_up = 0.5 * N_PILLARS * M_STAR * L_EDD_PER_MSUN     # upward half, erg/s
obs = load_pg0844()
nu_o = C / (obs['lam'] * ANGSTROM)
obs_nuLnu = nu_o * obs['host'] * NUL
i_v = int(np.argmin(np.abs(obs['lam'] - LAM_V)))
target = obs_nuLnu[i_v]
lam_o = obs['lam'][i_v]

print(f"{N_PILLARS} stars of {M_STAR:.0f} Msun at Eddington, upward half:")
print(f"  L_up = {L_up:.3e} erg/s")
print(f"observed constant component at {lam_o:.0f} A: nu L_nu = {target:.3e} erg/s")
print(f"  -> already {target/L_up:.1f}x the entire available luminosity\n")

# nu L_nu / L_bol for a blackbody at temperature T, evaluated at lam_o
nu = C / (lam_o * ANGSTROM)
T = np.logspace(2.5, 5.5, 4000)
frac = nu * b_nu(lam_o, T) / (SIGMA_SB * T**4 / np.pi)   # nu B_nu / (sigma T^4 / pi)
best = np.argmax(frac)
print(f"Blackbody fraction of L_bol emerging at {lam_o:.0f} A:")
print(f"  maximised at T = {T[best]:.0f} K, where nu L_nu / L_bol = {frac[best]:.3f}")
print(f"  best achievable nu L_nu = {frac[best]*L_up:.3e} erg/s")
print(f"  short of the observed component by a factor "
      f"{target/(frac[best]*L_up):.0f}\n")

# sphere radius corresponding to that optimal temperature, per pillar
L_one = 0.5 * M_STAR * L_EDD_PER_MSUN
R_opt = np.sqrt(L_one / (2.0 * np.pi * SIGMA_SB * T[best]**4))   # hemisphere
print(f"That optimal temperature corresponds to a hemisphere of radius "
      f"{R_opt/AU_CM:.3f} AU per pillar,")
print(f"i.e. shrinking the sphere all the way to {R_opt/AU_CM:.3f} AU still "
      f"falls {target/(frac[best]*L_up):.0f}x short.\n")

print("Total mass required instead, if every object radiates at Eddington and")
print("the emitter is tuned to the optimal temperature:")
m_req = target / frac[best] * 2.0 / L_EDD_PER_MSUN
print(f"  M_total = {m_req:.3e} Msun  ({m_req/N_PILLARS:.3e} Msun per pillar)")

print("\nHow the observed component varies with wavelength, against the ceiling:")
print(f"{'lam':>8}{'obs nuLnu':>13}{'max possible':>14}{'shortfall':>11}")
for lam, val in zip(obs['lam'], obs_nuLnu):
    nu_l = C / (lam * ANGSTROM)
    f = np.max(nu_l * b_nu(lam, T) / (SIGMA_SB * T**4 / np.pi))
    print(f"{lam:>8.0f}{val:>13.3e}{f*L_up:>14.3e}{val/(f*L_up):>10.0f}x")
