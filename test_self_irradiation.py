"""Estimate inner-disc self-irradiation vs lamp irradiation.

Keith's question: the viscous inner disc emits ionizing photons at a rate
locally comparable to the lamp (Fig A1). Is it safe to neglect the hot
flat inner disc irradiating the bowl and pillars further out?

Model: treat the inner disc as a Lambertian annular source with total
one-sided ionizing photon rate Q_disc = \int Phi_exact(r) 2 pi r dr.
Seen from a point at radius r and elevation z(r), the emission-direction
cosine is cos(theta_e) ~ (z - z_emit)/r (grazing), so the received flux is

    F_disc = Q_disc cos(theta_e) cos(theta_inc) / (pi r^2)

versus the isotropic lamp

    F_lamp = Q_lamp cos(theta_inc) / (4 pi d^2).

Ratio (same receiving surface, d ~ r):

    F_disc / F_lamp = 4 (Q_disc / Q_lamp) cos(theta_e).
"""
import numpy as np
from scipy.integrate import quad

# Constants (cgs)
H_PL = 6.626e-27
K_B = 1.381e-16
C_CGS = 2.998e10
NU_0 = 13.6 * 1.602e-12 / H_PL
G_CGS = 6.674e-8
M_SUN = 1.989e33
LD_CM = C_CGS * 86400.0

# Parameters (config_line.yaml)
M_BH = 7.0e7
R_ISCO_RG = 2.0
HLAMP = 0.05
LOG_PHI_INNER = 25.0
TV1 = 446.0
R0 = 20.0
H1 = 0.05
BETA = 10.0

r_g_ld = G_CGS * M_BH * M_SUN / C_CGS**2 / LD_CM
r_in = R_ISCO_RG * r_g_ld

# --- Lamp photon rate ---
d_in = np.sqrt(r_in**2 + HLAMP**2)
cos_in = HLAMP / d_in
Q_lamp = 10.0**LOG_PHI_INNER * 4.0 * np.pi * (d_in * LD_CM)**2 / cos_in
print(f'Q_lamp  = {Q_lamp:.3e} photons/s')

# --- Inner-disc photon rate: integrate exact Phi_visc over one side ---
def wien_integral(x0):
    if x0 > 100:
        return 0.0
    res, _ = quad(lambda u: u**2 / np.expm1(u), x0, max(x0 + 80.0, 100.0),
                  limit=200)
    return res

def phi_exact(r_ld):
    no_torque = np.sqrt(max(0.0, (r_ld - r_in) / r_ld))
    T = TV1 * (r_ld / R0)**(-0.75) * no_torque
    if T < 100:
        return 0.0
    x0 = H_PL * NU_0 / (K_B * T)
    pref = 2.0 * np.pi / C_CGS**2 * (K_B / H_PL)**3
    return pref * T**3 * wien_integral(x0)

# integrate Q_disc = int Phi(r) 2 pi r dr  (r in cm)
def integrand(r_ld):
    return phi_exact(r_ld) * 2.0 * np.pi * (r_ld * LD_CM) * LD_CM

Q_disc, _ = quad(integrand, r_in, 1.0, limit=400)  # emission beyond 1 ld negligible
print(f'Q_disc  = {Q_disc:.3e} photons/s (one side)')
print(f'Q_disc/Q_lamp = {Q_disc / Q_lamp:.3e}')

# --- Foreshortening cos(theta_e) at representative radii ---
def z_bowl(r_ld):
    return H1 * (r_ld / R0)**BETA

print('\n r [ld]   z(r) [ld]   cos(th_e)      F_disc/F_lamp')
for r in [2.0, 4.0, 10.0, 20.0]:
    z = z_bowl(r)
    cos_e = z / r                      # bowl surface
    ratio = 4.0 * (Q_disc / Q_lamp) * cos_e
    # pillar wall: top of pillar at z + h_p
    h_p = 0.04
    cos_e_pillar = (z + h_p) / r
    ratio_p = 4.0 * (Q_disc / Q_lamp) * cos_e_pillar
    print(f'  {r:5.1f}   {z:9.2e}   {cos_e:9.2e} (bowl)    {ratio:9.2e}')
    print(f'          {"":9}    {cos_e_pillar:9.2e} (pillar)  {ratio_p:9.2e}')
