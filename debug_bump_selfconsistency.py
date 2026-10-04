"""Is the pillar bump geometry self-consistent with the embedded-star luminosity?

In a radiation-hydrodynamic picture the bump is not a free shape: the embedded
star heats the gas, the gas puffs up to the new hydrostatic scale height, and
both the bump height and its lateral extent follow from L. This script checks
the fiducial geometry of test_fluxflux_embedded.py against that requirement:

  1. the ambient gas-pressure scale height H = c_s / Omega at the pillar belt,
  2. the temperature (hence L) needed to puff the disc to h_p = 0.5 ld,
  3. the heating radius R_heat where the star's flux equals the ambient disc
     flux, which is the natural lateral size of the heated region,
  4. the same comparison for a radiation-pressure supported bump.
"""
import numpy as np

from pillardisk.pillar_disk import C, DAY
import test_lag_spectrum as tls
from fluxflux_lib import DiscGrid, eta_lp, SIGMA_SB, L_EDD_PER_MSUN, AU_CM

G = 6.674e-8
M_SUN = 1.989e33
K_B = 1.381e-16
M_H = 1.673e-24
SIGMA_T_OVER_MP = 0.4                    # electron scattering opacity, cm^2/g
LD_CM = C * DAY

H_P, SIGMA_R, SIGMA_PHI = 0.5, 0.25, 0.30
HLAMP = 50 * 0.00399
R_BELT = 15.0                            # light days, middle of the pillar belt
M_BH = 7.0e7
M_EMB_PHYS = 7.37e3                      # Msun per pillar, total (from the scan)
MU = 0.6                                 # ionized mean molecular weight

# --- ambient disc temperature at the belt, from the actual model -----------
disk = tls.build_disk(hlamp=HLAMP)
disk.d = 287.0 * 3.086e24 / LD_CM
grid = DiscGrid(disk)
disk.tx_base *= (50.0 / eta_lp(grid)) ** 0.25
disk.tv_base *= 1.5
grid.tv2 *= 1.5
T2d = disk.get_temperature(grid.r2, grid.p2, compute_shadows=True)
ir = int(np.argmin(np.abs(disk.r - R_BELT)))
T_amb = float(T2d[ir].mean())
T_visc = float(np.interp(R_BELT, disk.r, disk.tv_base))
print(f"At r = {R_BELT} ld (bowl, no pillars): T_visc = {T_visc:.0f} K, "
      f"T_total = {T_amb:.0f} K")

# --- hydrostatic scale height ---------------------------------------------
r_cm = R_BELT * LD_CM
omega = np.sqrt(G * M_BH * M_SUN / r_cm**3)


def scale_height(T):
    c_s = np.sqrt(K_B * T / (MU * M_H))
    return c_s / omega


H_amb = scale_height(T_amb)
print(f"\nOmega = {omega:.3e} s^-1")
print(f"Ambient scale height H = {H_amb:.3e} cm = {H_amb/LD_CM:.4f} ld "
      f"= {H_amb/AU_CM:.1f} AU,  H/r = {H_amb/r_cm:.5f}")

h_p_cm = H_P * LD_CM
print(f"Model pillar height h_p = {H_P} ld = {h_p_cm:.3e} cm "
      f"= {h_p_cm/AU_CM:.0f} AU,  h_p / H_amb = {h_p_cm/H_amb:.1f}")

# temperature needed for gas pressure alone to support h_p
T_needed = T_amb * (h_p_cm / H_amb) ** 2
print(f"\nGas-pressure support of h_p needs T = {T_needed:.3e} K "
      f"({T_needed/T_amb:.0f}x the ambient temperature)")

# --- luminosity implied by that temperature over the bump footprint --------
for label, sp in [('sigma_phi=0.30', 0.30), ('sigma_phi=0.02', 0.02)]:
    a_bump = 2.0 * np.pi * (SIGMA_R * LD_CM) * (sp * r_cm)     # 2 pi sig_r sig_s
    l_needed = SIGMA_SB * T_needed**4 * a_bump
    print(f"  {label}: A_bump = {a_bump:.3e} cm^2 -> L needed = "
          f"{l_needed:.3e} erg/s = {l_needed/L_EDD_PER_MSUN:.3e} Msun at "
          f"Eddington")

# --- what the fiducial embedded luminosity actually does ------------------
L_emb = M_EMB_PHYS * L_EDD_PER_MSUN
print(f"\nFiducial embedded luminosity L = {L_emb:.3e} erg/s "
      f"({M_EMB_PHYS:.3g} Msun at Eddington)")
R_heat = np.sqrt(L_emb / (4.0 * np.pi * SIGMA_SB * T_amb**4))
print(f"Heating radius where L/(4 pi R^2) = sigma T_amb^4: "
      f"R_heat = {R_heat:.3e} cm = {R_heat/AU_CM:.1f} AU "
      f"= {R_heat/LD_CM:.4f} ld")
T_heat_mean = (L_emb / (SIGMA_SB * np.pi * R_heat**2)) ** 0.25
print(f"Mean temperature over that heated area: {T_heat_mean:.0f} K "
      f"= {T_heat_mean/T_amb:.2f} x T_amb")
H_heat = scale_height(T_heat_mean)
print(f"Scale height at that temperature: {H_heat/LD_CM:.4f} ld "
      f"= {H_heat/AU_CM:.1f} AU (vs h_p = {H_P} ld)")

# --- radiation-pressure puffing -------------------------------------------
# vertical radiation force per gram kappa F / c vs vertical gravity Omega^2 z
print("\nRadiation pressure check (electron scattering, kappa = 0.4 cm^2/g):")
for R_au in [70.0, 284.0]:
    R = R_au * AU_CM
    F_rad = L_emb / (4.0 * np.pi * R**2)
    a_rad = SIGMA_T_OVER_MP * F_rad / C
    z_eq = a_rad / omega**2                     # height where gravity balances
    print(f"  at R = {R_au:.0f} AU: F = {F_rad:.3e} erg/cm^2/s, "
          f"a_rad/g at z=H_amb is {a_rad/(omega**2*H_amb):.2f}, "
          f"balance height z = {z_eq/LD_CM:.4f} ld")
