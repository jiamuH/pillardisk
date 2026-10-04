"""Referee reply, Section 7 bullet 8: where should the two pillar hotspots sit in
the velocity-delay plane, and where do they sit in the computed map?

Geometry, as coded in pillar_line_time_cloudy.compute_velocity_delay_map_cloudy:
    v_los = + v_K(r) sin(i) sin(phi)          (positive = receding = redshift)
    tau   = d_lamp - (r cos(phi) sin(i) + h cos(i)),
    d_lamp = sqrt(r^2 + (h - h_lamp)^2),   distances in light days = days.

Pillars: r_p = 4 light days, phi_p = 3 pi / 4 and 0, height 0.15 light days,
inclination 45 degrees, M_BH = 7e7 Msun, lamp height 0.05 light days.

Run:  python3 debug_hotspot_positions.py
"""
import numpy as np

G = 6.674e-8
M_SUN = 1.989e33
C = 2.998e10
LIGHT_DAY = C * 86400.
LAMBDA0 = 6562.8

M_BH = 7e7 * M_SUN
R_P, H_P, H_LAMP = 4.0, 0.15, 0.05
H0, R0, BETA = 0.05, 20.0, 10.0
INCL = np.radians(45.)

v_k = np.sqrt(G * M_BH / (R_P * LIGHT_DAY)) / 1e5          # km/s
period = 2. * np.pi * R_P * LIGHT_DAY / (v_k * 1e5) / 86400.  # days
print('Keplerian speed at r_p = 4 light days: %.0f km/s' % v_k)
print('Orbital period: %.0f days' % period)
print('Maximum projected speed v_K sin i: %.0f km/s\n' % (v_k * np.sin(INCL)))

for label, phi in (('phi_p = 3 pi / 4', 3. * np.pi / 4.), ('phi_p = 0', 0.)):
    v_los = v_k * np.sin(INCL) * np.sin(phi)
    print(label)
    print('  line-of-sight velocity %+.0f km/s  (%s), wavelength %.1f A'
          % (v_los, 'redshift' if v_los > 1 else 'line centre' if abs(v_los) < 1 else 'blueshift',
             LAMBDA0 * (1. + v_los * 1e5 / C)))
    for where, h in (('disc surface at the pillar', H0 * (R_P / R0) ** BETA),
                     ('pillar top', H0 * (R_P / R0) ** BETA + H_P)):
        d_lamp = np.sqrt(R_P ** 2 + (h - H_LAMP) ** 2)
        tau = d_lamp - (R_P * np.cos(phi) * np.sin(INCL) + h * np.cos(INCL))
        print('  delay at the %-28s %.2f days' % (where + ':', tau))

# Where are the hotspots in the computed map? Use the excess of the
# with-pillars map over the no-pillar map.
d = np.load('plots/vdm_manual_pillars_maps.npz')
lam, tau, psi, psi_no = d['lam'], d['tau'], d['psi'], d['psi_no']
excess = psi - psi_no
print('\nStrongest excess response (with pillars minus without) in the computed map:')
for label, lam_lo, lam_hi in (('red side, 6600 to 6800 A', 6600., 6800.),
                              ('around line centre, 6520 to 6600 A', 6520., 6600.),
                              ('blue side, 6300 to 6520 A', 6300., 6520.)):
    sel = (lam >= lam_lo) & (lam <= lam_hi)
    sub = excess[:, sel]
    it, il = np.unravel_index(np.argmax(sub), sub.shape)
    lam_pk = lam[sel][il]
    print('  %-36s peak at %.1f A = %+.0f km/s, tau = %.1f days, excess = %.2e'
          % (label + ':', lam_pk, (lam_pk / LAMBDA0 - 1.) * C / 1e5, tau[it], sub[it, il]))
print('Delay axis of the map: 0 to %.0f days in %d bins' % (tau[-1], len(tau)))
