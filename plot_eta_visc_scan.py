"""Mechanism A: the constant viscous floor as a fake host, vs lamp strength.

The viscous emission never varies (viscous timescale >> a monitoring campaign),
so a flux-flux analysis dumps it into the constant ("host") component. Here we
score the recovered host against the STELLAR host alone (HOST, not HOST+floor),
so that contamination shows up, and sweep ETA_LP downward (weaker lamp -> bigger
viscous fraction). Pillar (with shadows) vs bowl (no shadows) tests whether
shadows feed it. Fixed: lamp 50 r_g, thin pillars, 20% campaign, current host
(fraction 0.2 at V, beta=3).

Two curves per disc at V band:
  (HOST+floor)/HOST   = ideal flux-flux overestimate = pure viscous contamination
  recovered/HOST      = what the anchor method actually returns (floor + curvature)
"""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import C, DAY, H, K, ANGSTROM
from test_lag_spectrum import build_disk, add_pillars

H_P, SIGMA_R, SIGMA_PHI = 0.5, 0.25, 0.3
HLAMP = 50 * 0.00399
HOST_BETA, HOST_5100 = 3.0, 1.825
bands = np.array([1928., 2246., 2600., 3465., 4392., 5468., 6215., 7545., 8700.])
iV = int(np.argmin(np.abs(bands - 5100.0)))
HOST = HOST_5100 * (bands / 5100.0) ** HOST_BETA
s_vals = np.linspace(1.08, 1.35, 12)               # 20% campaign window
eta_vals = np.array([1.0, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0])


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


def metrics(with_pillars, eta):
    dk = build_disk(hlamp=HLAMP)
    if with_pillars:
        add_pillars(dk, 100, H_P, SIGMA_R, SIGMA_PHI, r_min=10.0)
    dk.tx_base *= (eta / eta_lp(dk))**0.25
    F, floor = band_fluxes(dk, s_vals)
    FS = 9.0 / F[:, iV].max()
    F *= FS; floor *= FS
    rec = recovered(F + HOST)
    ideal = (HOST + floor) / HOST                  # pure viscous contamination
    actual = rec / HOST                            # anchor recovery vs stellar host
    return ideal, actual, floor


ideal_p = np.zeros((len(eta_vals), len(bands)))
ideal_b = np.zeros((len(eta_vals), len(bands)))
act_p = np.zeros((len(eta_vals), len(bands)))
act_b = np.zeros((len(eta_vals), len(bands)))
print(f"{'eta':>5} {'floorV/host':>11} {'(h+fl)/h V':>10} {'rec/h V':>8}  (pillar | bowl)")
for i, eta in enumerate(eta_vals):
    ideal_p[i], act_p[i], fp = metrics(True, eta)
    ideal_b[i], act_b[i], fb = metrics(False, eta)
    print(f"{eta:5.0f} {fp[iV]/HOST[iV]:5.2f}/{fb[iV]/HOST[iV]:<5.2f} "
          f"{ideal_p[i, iV]:4.2f}/{ideal_b[i, iV]:<5.2f} "
          f"{act_p[i, iV]:4.2f}/{act_b[i, iV]:<5.2f}")

# --- plot (V band) ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(eta_vals, ideal_p[:, iV], '-o', color='purple', lw=3, ms=7,
        label=r'$\rm pillar:~(host+floor)/host~(viscous~only)$')
ax.plot(eta_vals, ideal_b[:, iV], '--s', color='purple', lw=2.2, ms=6, alpha=0.6,
        label=r'$\rm bowl:~(host+floor)/host$')
ax.plot(eta_vals, act_p[:, iV], '-o', color='orangered', lw=3, ms=7,
        label=r'$\rm pillar:~recovered/host~(anchor)$')
ax.plot(eta_vals, act_b[:, iV], '--s', color='orangered', lw=2.2, ms=6, alpha=0.6,
        label=r'$\rm bowl:~recovered/host~(anchor)$')
ax.axhline(1.0, color='gray', ls=':', lw=2, alpha=0.8)
ax.set_xscale('log')
ax.set_xlabel(r'$\eta_{\rm LP}~\rm (lamp/viscous~power)$', fontsize=18)
ax.set_ylabel(r'$\rm recovered~host~/~stellar~host~(V)$', fontsize=18)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=11, frameon=False, loc='best')
plt.tight_layout()
plt.savefig('plots/test_eta_visc_scan.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_eta_visc_scan.png")
