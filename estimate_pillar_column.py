"""Referee reply, Section 4 bullet 7: how large is the hydrogen column of a
pillar, that is of the disc spot heated by an embedded star?

The paper states no disc gas density, so the disc is taken to be marginally
self-gravitating (Toomre Q = 1), as in the self-gravitating outer-disc models
the paper cites. Then

    Sigma = c_s * Omega / (pi * G),     H = c_s / Omega,

with the sound speed c_s from the disc temperature. Columns quoted are
hydrogen columns, N_H = X * Sigma / m_H with hydrogen mass fraction X = 0.7.

Geometry and temperature follow the paper: M_BH = 7e7 Msun, viscous
temperature T(r) = 446 K * (r / 20 light days)^-0.75, pillars of height
h_p = 0.04 light days and radial width sigma_r = 0.1 light days at
r = 10 light days (the N = 100 population), and h_p = 0.15, sigma_r = 0.5 at
r = 4 light days (the single-pillar case).

Run:  python3 estimate_pillar_column.py
"""
import numpy as np

G = 6.674e-8          # cm^3 g^-1 s^-2
K_B = 1.381e-16       # erg K^-1
M_H = 1.673e-24       # g
M_SUN = 1.989e33      # g
LIGHT_DAY = 2.59e15   # cm
X_HYDROGEN = 0.7
MU = 1.3              # mean molecular weight, neutral atomic gas with helium
SIGMA_912 = 6.3e-18   # cm^2, hydrogen photoionization cross-section at the edge

M_BH = 7e7 * M_SUN
T_REFERENCE, R_REFERENCE, ALPHA = 446., 20. * LIGHT_DAY, 0.75

CASES = [('single pillar', 4., 0.15, 0.5),
         ('N = 100 pillar', 10., 0.04, 0.1)]


def disc_structure(r_cm):
    T = T_REFERENCE * (r_cm / R_REFERENCE) ** (-ALPHA)
    c_s = np.sqrt(K_B * T / (MU * M_H))
    omega = np.sqrt(G * M_BH / r_cm ** 3)
    sigma_gas = c_s * omega / (np.pi * G)     # Toomre Q = 1
    scale_height = c_s / omega
    density = sigma_gas / (2. * scale_height)  # g cm^-3, midplane
    n_h = X_HYDROGEN * density / M_H           # hydrogen nuclei cm^-3
    return T, c_s, omega, sigma_gas, scale_height, n_h


print('Marginally self-gravitating disc (Toomre Q = 1), M_BH = 7e7 Msun\n')
for label, r_ld, h_p_ld, sigma_r_ld in CASES:
    r_cm = r_ld * LIGHT_DAY
    T, c_s, omega, sigma_gas, H, n_h = disc_structure(r_cm)
    vertical = X_HYDROGEN * sigma_gas / M_H
    through_height = n_h * h_p_ld * LIGHT_DAY
    through_width = n_h * sigma_r_ld * LIGHT_DAY
    print('%s at r = %.0f light days' % (label, r_ld))
    print('  disc temperature                     %.0f K' % T)
    print('  sound speed                          %.2e cm/s' % c_s)
    print('  surface density Sigma                %.2e g cm^-2' % sigma_gas)
    print('  scale height H                       %.2e cm  (%.2e light days)'
          % (H, H / LIGHT_DAY))
    print('  pillar height h_p                    %.2e cm  (%.2f light days, %.0f H)'
          % (h_p_ld * LIGHT_DAY, h_p_ld, h_p_ld * LIGHT_DAY / H))
    print('  midplane hydrogen density            %.2e cm^-3' % n_h)
    print('  vertical column through the disc     log N_H = %.1f' % np.log10(vertical))
    print('  column over the pillar height        log N_H = %.1f (at midplane density)'
          % np.log10(through_height))
    print('  column across the pillar width       log N_H = %.1f (at midplane density)'
          % np.log10(through_width))
    print()

print('Geometric column across a pillar at the densities of the Cloudy grid')
print('(path length = radial width sigma_r, no assumption about the disc):')
print('%12s %22s %22s' % ('log n_H', 'log N_H, sigma_r = 0.1 ld', 'log N_H, sigma_r = 0.5 ld'))
for log_n in (9., 10., 11., 12.):
    n = 10. ** log_n
    print('%12.1f %22.1f %22.1f'
          % (log_n, np.log10(n * 0.1 * LIGHT_DAY), np.log10(n * 0.5 * LIGHT_DAY)))
print()

print('For comparison:')
print('  Cloudy stopping column in the grid      log N_H = 25.0')
print('  column where an ionizing photon at 912 A reaches optical depth 1,')
print('  if the gas were neutral                 log N_HI = %.1f'
      % np.log10(1. / SIGMA_912))
