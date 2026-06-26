"""Does thinning the pillars reduce the optical/UV bump? Scan the width.

Thinning the pillars shrinks the illuminated-wall AREA but steepens the
lamp-facing wall (hotter walls), so the net effect on the bump is not obvious.
Here we scale sigma_r and sigma_phi together by a factor, rebuild the pillar
disc at each width (fixed pillar positions via the default seed), and measure
the bump as the bright-state pillar/bowl flux ratio at a few wavelengths. A
ratio above 1 is the bump; ratio -> 1 means the bump has gone away.

H_P is lowered to 0.5 for this test. Each disc is renormalized to the same
bolometric eta_LP = 3, so the comparison is at fixed total irradiation budget
and isolates how the wall geometry redistributes that budget into the blue.
"""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import C, DAY, H, K, ANGSTROM
from test_lag_spectrum import build_disk, add_pillars

HLAMP = 10 * 0.00399          # lamp height = 10 r_g (match the mock)
H_P = 0.5                     # lowered for this bump test
ETA_LP = 3.0                  # current irradiation strength
S_BRIGHT = 1.35               # bright lamp state for the SED
SR0, SP0 = 0.5, 0.6           # current widths (scaled together by f)
wl_probe = np.array([2000.0, 3000.0, 5000.0])   # UV, near-UV, optical


def b_nu(lam_A, T):
    nu = C / (lam_A * ANGSTROM)
    x = np.clip(H * nu / (K * np.maximum(T, 1.0)), 1e-10, 700.0)
    return 2.0 * H * nu**3 / C**2 / (np.exp(x) - 1.0)


def eta_lp(dk):
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    ds = np.sqrt(np.diff(dk.r)[:, None]**2 + (h2[1:, :] - h2[:-1, :])**2)
    da = ds * dk.r[:-1, None] * dk.dphi
    tv2 = np.interp(r2.ravel(), dk.r, dk.tv_base).reshape(r2.shape)
    T = dk.get_temperature(r2, p2, compute_shadows=True)
    tx4 = np.maximum(T**4 - tv2**4, 0.0)
    return np.sum(tx4[:-1, :] * da) / np.sum(tv2[:-1, :]**4 * da)


def disc_sed(dk, wl, s):
    """Relative F_nu (arbitrary scale) for disc dk at lamp state s; the absolute
    scale cancels in the pillar/bowl ratio."""
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    dsg = np.sqrt(np.diff(dk.r)[:, None]**2 + (h2[1:, :] - h2[:-1, :])**2)
    px = -(h2[1:, :] - h2[:-1, :]) / (dsg + 1e-10) * np.cos(dk.phi)[None, :]
    pz = np.diff(dk.r)[:, None] / (dsg + 1e-10)
    dot = np.where(dk.ex * px + dk.ez * pz > 0, dk.ex * px + dk.ez * pz, 0.0)
    wg = dsg * dk.r[:-1, None] * dk.dphi * dot
    t0 = dk.tx_base.copy(); dk.tx_base = t0 * s
    T = dk.get_temperature(r2, p2, compute_shadows=True); dk.tx_base = t0
    return np.array([np.sum(b_nu(l, T[:-1, :] * dk.fcol) / dk.fcol**4 * wg)
                     for l in wl])


# bowl reference (no pillars), renormalized to eta_LP = 3
bowl = build_disk(hlamp=HLAMP)
bowl.tx_base *= (ETA_LP / eta_lp(bowl))**0.25
bowl_sed = disc_sed(bowl, wl_probe, S_BRIGHT)

f_vals = np.array([0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4])
sr_vals = SR0 * f_vals
ratios = np.zeros((len(f_vals), len(wl_probe)))
print(f"{'sig_r':>6} {'sig_phi':>7} " + " ".join(f"{w:>5.0f}A" for w in wl_probe))
for i, f in enumerate(f_vals):
    dk = build_disk(hlamp=HLAMP)
    add_pillars(dk, 100, H_P, SR0 * f, SP0 * f, r_min=10.0)
    dk.tx_base *= (ETA_LP / eta_lp(dk))**0.25
    ratios[i] = disc_sed(dk, wl_probe, S_BRIGHT) / bowl_sed
    print(f"{SR0*f:6.3f} {SP0*f:7.3f} " + " ".join(f"{r:6.3f}" for r in ratios[i]))

# --- plot ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

clr = {2000: 'royalblue', 3000: 'seagreen', 5000: 'orangered'}
fig, ax = plt.subplots(figsize=(9, 7))
for j, w in enumerate(wl_probe):
    ax.plot(sr_vals, ratios[:, j], '-o', color=clr[int(w)], lw=3, ms=7,
            label=r'$' + f'{int(w)}' + r'~\AA$')
ax.axhline(1.0, color='gray', ls=':', lw=2, alpha=0.8)         # no bump
ax.axvline(SR0, color='k', ls=':', lw=1.5, alpha=0.6)          # current sigma_r
ax.text(SR0 + 0.01, ax.get_ylim()[0], r'$\rm current$', rotation=90,
        va='bottom', ha='left', fontsize=12)
ax.set_xlabel(r'$\sigma_r~\rm [light~days]~(\sigma_\phi=1.2\,\sigma_r)$', fontsize=18)
ax.set_ylabel(r'$\nu L_\nu~\rm (pillar)~/~\nu L_\nu~\rm (bowl)$', fontsize=18)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=14, frameon=False, title=r'$\rm wavelength$')
plt.tight_layout()
plt.savefig('plots/test_bump_width_scan.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_bump_width_scan.png")
