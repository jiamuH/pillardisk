"""Host-recovery bias vs lamp-post height (old anchor method).

Vary the lamp height, holding eta_LP = 3 fixed (so the bolometric irradiation
budget is constant and we isolate the geometric effect of moving the lamp up).
Current pillars (H_P=0.5, sigma_r=0.25, sigma_phi=0.3) and the 20% campaign
window. For each height: sweep the lamp, add the host, run the old recovery
(driver -> per-band line -> anchor at the bluest band's F=0), and compare the
recovered host to the true host (input host + viscous floor). Pillar = solid,
bowl = dashed.
"""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import C, DAY, H, K, ANGSTROM
from test_lag_spectrum import build_disk, add_pillars

H_P, SIGMA_R, SIGMA_PHI = 0.5, 0.25, 0.3       # current thin pillars
ETA_LP = 3.0
R_G = 0.00399                                   # light days, for labels
bands = np.array([1928., 2246., 2600., 3465., 4392., 5468., 6215., 7545., 8700.])
iV = int(np.argmin(np.abs(bands - 5100.0)))
HOST_BETA, HOST_5100 = 2.0, 3.5
HOST = HOST_5100 * (bands / 5100.0) ** HOST_BETA
s_vals = np.linspace(1.08, 1.35, 12)            # 20% campaign window
hlamp_ld = np.array([0.04, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0])


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


def band_fluxes(dk, s_vals):
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    dh = h2[1:, :] - h2[:-1, :]
    ds = np.sqrt(np.diff(dk.r)[:, None]**2 + dh**2)
    px = -(dh / (ds + 1e-10)) * np.cos(dk.phi)[None, :]
    pz = np.diff(dk.r)[:, None] / (ds + 1e-10)
    dot = np.where(dk.ex * px + dk.ez * pz > 0, dk.ex * px + dk.ez * pz, 0.0)
    wg = ds * dk.r[:-1, None] * dk.dphi * dot
    tv2 = np.interp(r2.ravel(), dk.r, dk.tv_base).reshape(r2.shape)
    t0 = dk.tx_base.copy()
    F = np.zeros((len(s_vals), len(bands)))
    for j, s in enumerate(s_vals):
        dk.tx_base = t0 * s
        T = dk.get_temperature(r2, p2, compute_shadows=True)
        F[j] = np.array([np.sum(b_nu(l, T[:-1, :] * dk.fcol) / dk.fcol**4 * wg)
                         for l in bands])
    dk.tx_base = t0
    floor = np.array([np.sum(b_nu(l, tv2[:-1, :] * dk.fcol) / dk.fcol**4 * wg)
                      for l in bands])
    return F, floor


def driving_lc(F):
    B = F.mean(axis=0); A = F.std(axis=0)
    X = ((F - B) / A).mean(axis=1)
    return (X - X.mean()) / X.std()


def recovered(F_h):
    X = driving_lc(F_h)
    fit = np.array([np.polyfit(X, F_h[:, k], 1) for k in range(len(bands))])
    anchor = np.max(-fit[:, 1] / fit[:, 0])
    return fit[:, 0] * anchor + fit[:, 1]


def host_ratio(hl, with_pillars):
    dk = build_disk(hlamp=hl)
    if with_pillars:
        add_pillars(dk, 100, H_P, SIGMA_R, SIGMA_PHI, r_min=10.0)
    dk.tx_base *= (ETA_LP / eta_lp(dk))**0.25
    F, floor = band_fluxes(dk, s_vals)
    FS = 9.0 / F[:, iV].max()
    F *= FS; floor *= FS
    rec = recovered(F + HOST)
    return rec / (HOST + floor)


show = [int(np.argmin(np.abs(bands - w))) for w in (3465, 5468, 8700)]
rp = np.zeros((len(hlamp_ld), len(bands)))
rb = np.zeros((len(hlamp_ld), len(bands)))
print(f"{'hlamp[ld]':>9} {'[r_g]':>6}  pillar(3465/5468/8700) | bowl")
for i, hl in enumerate(hlamp_ld):
    rp[i] = host_ratio(hl, True)
    rb[i] = host_ratio(hl, False)
    pp = "/".join(f"{rp[i, k]:.2f}" for k in show)
    bb = "/".join(f"{rb[i, k]:.2f}" for k in show)
    print(f"{hl:9.2f} {hl/R_G:6.0f}  {pp} | {bb}")

# --- plot ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

clr = {3465: 'royalblue', 5468: 'seagreen', 8700: 'orangered'}
fig, ax = plt.subplots(figsize=(9, 7))
for k, w in zip(show, (3465, 5468, 8700)):
    ax.plot(hlamp_ld, rp[:, k], '-o', color=clr[w], lw=3, ms=7,
            label=r'$\rm pillar,~' + f'{w}' + r'~\AA$')
    ax.plot(hlamp_ld, rb[:, k], '--s', color=clr[w], lw=2.2, ms=6, alpha=0.7,
            label=r'$\rm bowl,~' + f'{w}' + r'~\AA$')
ax.axhline(1.0, color='gray', ls=':', lw=2, alpha=0.8)         # perfect / overestimate above
ax.axvline(0.04, color='k', ls=':', lw=1.5, alpha=0.6)         # current
ax.text(0.043, ax.get_ylim()[0], r'$\rm current$', rotation=90,
        va='bottom', ha='left', fontsize=12)
ax.set_xscale('log')
ax.set_xlabel(r'$\rm lamp~height~[light~days]$', fontsize=18)
ax.set_ylabel(r'$f_\nu^{\rm recovered}/f_\nu^{\rm true}~\rm (host)$', fontsize=18)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=11, frameon=False, ncol=2, loc='best')
plt.tight_layout()
plt.savefig('plots/test_lampheight_scan.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_lampheight_scan.png")
