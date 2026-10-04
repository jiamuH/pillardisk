"""One embedded star per pillar: can a few hundred stars supply the observed
constant ("host") component of PG 0844+349?

Earlier runs treated the mass per pillar as a summed embedded population. If
each pillar is a single star, that number is the stellar mass, so the budget has
to be redone. This integrates the observed constant component to get the
luminosity that has to be supplied, then asks what one star per pillar can
actually deliver.
"""
import numpy as np

from fluxflux_lib import load_pg0844, L_EDD_PER_MSUN
from pillardisk.pillar_disk import C, ANGSTROM

MPC_CM = 3.086e24
DL_CM = 287.0 * MPC_CM
L_SUN = 3.828e33

obs = load_pg0844()
lam, host = obs['lam'], obs['host']
nu = C / (lam * ANGSTROM)

# nu L_nu and the integrated luminosity over the observed range. F_nu is in mJy
# (1e-26 erg s^-1 cm^-2 Hz^-1); integrate in frequency, which the descending
# wavelength order makes increasing in nu.
NUL = 1e-26 * 4.0 * np.pi * DL_CM**2
order = np.argsort(nu)
L_int = np.trapz(host[order] * 1e-26 * 4.0 * np.pi * DL_CM**2, nu[order])

print("Observed constant ('host') component of PG 0844+349:")
for l, h in zip(lam, host):
    print(f"  {l:7.0f} A : F_nu = {h:6.3f} mJy, nu L_nu = "
          f"{C/(l*ANGSTROM)*h*NUL:.3e} erg/s")
print(f"\nIntegrated over the observed range "
      f"({lam.min():.0f}-{lam.max():.0f} A): L = {L_int:.3e} erg/s")
print("  (a lower bound: the component is still rising redward of 8177 A)")
# crude bolometric correction for a blackbody peaking near the red end
for bc in [1.5, 2.0, 3.0]:
    print(f"  x{bc:.1f} bolometric correction -> {L_int*bc:.3e} erg/s")

L_NEED = 2.0 * L_int      # adopt a factor 2 bolometric correction
print(f"\nAdopting L_need = {L_NEED:.3e} erg/s.")
print("Only the upward half of a buried star's luminosity is observable, so the")
print("population must supply 2 x L_need in total Eddington luminosity.")
L_SUPPLY_NEED = 2.0 * L_NEED

print("\nWhat one star per pillar delivers (L_Edd each, half escaping upward):")
print(f"{'N stars':>9}{'m_star':>10}{'M_total':>12}{'L_Edd,tot':>12}"
      f"{'L_up':>12}{'fraction of need':>18}")
for n in [100, 300, 500]:
    for m in [100.0, 200.0, 500.0]:
        l_tot = n * m * L_EDD_PER_MSUN
        print(f"{n:>9}{m:>10.0f}{n*m:>12,.0f}{l_tot:>12.3e}"
              f"{0.5*l_tot:>12.3e}{0.5*l_tot/L_NEED*100:>17.2f}%")

print("\nStellar mass required per pillar to supply the whole component:")
print(f"{'N stars':>9}{'m_star needed':>16}{'M_total':>14}")
for n in [100, 300, 500]:
    m_req = L_SUPPLY_NEED / (n * L_EDD_PER_MSUN)
    print(f"{n:>9}{m_req:>16,.0f}{n*m_req:>14,.0f}")

print("\nFor reference, the most massive stars:")
for m in [100.0, 200.0, 500.0]:
    print(f"  m = {m:.0f} Msun: L_Edd = {m*L_EDD_PER_MSUN:.3e} erg/s "
          f"= {m*L_EDD_PER_MSUN/L_SUN:.3e} Lsun")
