"""Reproduce the Cai+24 Flux Variation Gradient construction for the BOWL.

Their method (their Figure 2, right panels): in the two-band flux-flux plane
(g-band flux on x, u-band flux on y), the variable nucleus traces a straight
line. Fit that line to the observed epochs. The steady host galaxy is the point
where this nucleus-variation line crosses a line through the origin whose slope
is the host galaxy's own u/g color (from the assumed host spectrum). That
crossing is the recovered host; compare it to the true (input) host.

The bowl has nearly constant nucleus color (little curvature), so this should
recover the host well, with only the small overestimate the literature reports.

Imports the already-computed bowl arrays from test_fluxflux_mock.
"""
import numpy as np
import matplotlib.pyplot as plt
import test_fluxflux_mock as m

# --- pick the two bands: u-like (blue) and g-like (red) ---
i_u = int(np.argmin(np.abs(m.bands - 3465.0)))     # blue band  (SDSS u analog)
i_g = int(np.argmin(np.abs(m.bands - 4392.0)))     # red  band  (SDSS g analog)
lam_u, lam_g = m.bands[i_u], m.bands[i_g]

# --- bowl fluxes (mJy), already scaled. total = nucleus + host ---
tot = m.F_bowl_h                # AGN + host, all lamp states
agn = m.F_bowl                  # pure AGN (no host)
host_u, host_g = m.HOST[i_u], m.HOST[i_g]
camp = m.camp                   # observed campaign states

xg_tot, yu_tot = tot[:, i_g], tot[:, i_u]
xg_agn, yu_agn = agn[:, i_g], agn[:, i_u]

# --- fit the nucleus-variation line to the OBSERVED total-flux epochs ---
m_tot, c_tot = np.polyfit(xg_tot[camp], yu_tot[camp], 1)   # y = m_tot x + c_tot
m_agn, c_agn = np.polyfit(xg_agn[camp], yu_agn[camp], 1)   # pure-AGN line

# --- host-color line through the origin: slope = host u/g flux ratio ---
slope_host = host_u / host_g                                # y = slope_host x

# --- recovered host = crossing of the nucleus line and the host-color line ---
xg_rec = c_tot / (slope_host - m_tot)
yu_rec = slope_host * xg_rec
print(f"bands: u={lam_u:.0f} A, g={lam_g:.0f} A")
print(f"nucleus-variation slope (u/g) = {m_tot:.3f}, host-color slope = {slope_host:.3f}")
print(f"recovered host:  g={xg_rec:.3f}  u={yu_rec:.3f} mJy")
print(f"true     host:  g={host_g:.3f}  u={host_u:.3f} mJy")
print(f"overestimate factor:  g={xg_rec/host_g:.2f}x  u={yu_rec/host_u:.2f}x")

# --- plot ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

fig, ax = plt.subplots(figsize=(8, 7.5))

# total-flux epochs: observed campaign opaque, unobserved faint tail transparent
ax.plot(xg_tot, yu_tot, 'o', color='k', ms=6, alpha=0.18)
ax.plot(xg_tot[camp], yu_tot[camp], 'o', color='k', ms=7, alpha=0.9,
        label=r'$\rm total~(nucleus+host),~observed$')
# pure-nucleus epochs (pass through the origin)
ax.plot(xg_agn[camp], yu_agn[camp], 'x', color='royalblue', ms=8, mew=2.0,
        label=r'$\rm pure~nucleus$')

# line span for drawing
xmax = 1.05 * xg_tot.max()
xline = np.array([0.0, xmax])
# nucleus-variation line (fit to observed total flux), extended to the crossing
ax.plot(xline, m_tot * xline + c_tot, '-', color='k', lw=2.0,
        label=r'$\rm nucleus~variation~line~(fit)$')
# pure-nucleus line (through origin), light blue dashed
ax.plot(xline, m_agn * xline + c_agn, '--', color='royalblue', lw=1.6, alpha=0.8)
# host-color line through the origin
ax.plot(xline, slope_host * xline, '--', color='firebrick', lw=2.0,
        label=r'$\rm host~color~line~(through~origin)$')

# recovered host = crossing point
ax.plot(xg_rec, yu_rec, 'o', mfc='white', mec='firebrick', mew=2.5, ms=14,
        label=r'$\rm recovered~host~(crossing)$')
# true (input) host
ax.plot(host_g, host_u, 'D', color='royalblue', ms=12, mec='k', mew=0.8,
        label=r'$\rm true~host~(input)$')

ax.set_xlim(0, xmax)
ax.set_ylim(0, 1.05 * yu_tot.max())
ax.set_xlabel(r'$f_\nu~\rm [mJy]~(g\!-\!band,~' + f'{lam_g:.0f}' + r'~\AA)$',
              fontsize=18)
ax.set_ylabel(r'$f_\nu~\rm [mJy]~(u\!-\!band,~' + f'{lam_u:.0f}' + r'~\AA)$',
              fontsize=18)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=12, frameon=False, loc='upper left')
plt.tight_layout()
plt.savefig('plots/test_fvg_bowl.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved plots/test_fvg_bowl.png")
