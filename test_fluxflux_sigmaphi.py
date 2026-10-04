"""Pillar azimuthal width scan for the embedded-star hot spot picture.

sigma_phi is an angle, so at r ~ 15 light days the fiducial sigma_phi = 0.3 rad
smears each pillar into a 4.5 ld arc against its 0.25 ld radial width, an 18:1
aspect ratio. This script repeats the embedded-star flux-flux experiment of
test_fluxflux_embedded.py for a range of sigma_phi, from that arc down to a
nearly round spot, holding the embedded object mass per pillar fixed.

Because the spot luminosity is fixed while its footprint area A_p shrinks with
sigma_phi, the spot temperature rises as T_spot ~ (L/A_p)^(1/4): the constant
component keeps its luminosity but moves blueward and becomes more compact.

Outputs, in plots/ :
  test_fluxflux_sigmaphi_tmap_s<value>.png   line-of-sight temperature map
  test_fluxflux_sigmaphi_host.png            recovered "host" for every
                                             sigma_phi, with the input spot
                                             blackbodies and the PG 0844 data

The azimuthal grid is refined to NPHI so the narrowest pillar is resolved:
sigma_phi = 0.02 needs dphi well below that, i.e. nphi of order 1440.
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pillardisk.pillar_disk import C, ANGSTROM, PC_TO_LD
import test_lag_spectrum as tls
from fluxflux_lib import (LD_CM, DiscGrid, spot_t4, diffuse_amplitude,
                          diffuse_flux, lamp_states, eta_lp, driving_lc,
                          nonvar_anchor, load_pg0844, overlay_pg0844, PG_KEYS,
                          mass_label, ticks, los_map, L_EDD_PER_MSUN, AU_CM)

# --- geometry, identical to test_fluxflux_embedded.py except sigma_phi ------
H_P, SIGMA_R = 0.5, 0.25
SIGMA_PHI_SCAN = [0.30, 0.10, 0.05, 0.02]      # radians
SIGMA_PHI_REF = 0.30                            # anchors the flux normalization
HLAMP = 50 * 0.00399
N_PILLARS = 100
DL_MPC = 287.0
ETA_LP = 50.0
NPHI = 1440                    # refined so sigma_phi = 0.02 spans ~4.6 cells

# `--quick` runs a coarse two-width version for debugging, not for results
if '--quick' in sys.argv:
    SIGMA_PHI_SCAN, NPHI = [0.30, 0.10], 360
    print("QUICK MODE: coarse grid, two widths only, results are not final")

# --- embedded stars: TOTAL mass per pillar, summed over the population ------
M_EMBEDDED = 5.0e4             # Msun per pillar, in model (pre-FLUX_SCALE) units
M_STAR_UNIT = 200.0

# --- bands, lamp sweep, diffuse glow and true host, as in the parent script --
bands = np.array([1928.0, 2246.0, 2600.0, 3465.0, 4392.0,
                  5468.0, 6215.0, 7545.0, 8700.0])
iV = int(np.argmin(np.abs(bands - 5100.0)))
F_DIFF = 0.3
DNU = C / (bands.min() * ANGSTROM) - C / (bands.max() * ANGSTROM)
S_FULL, CAMP_FRAC = 1.35, 0.20
S_LO = S_FULL * (1.0 - CAMP_FRAC)
s_vals = np.concatenate([np.linspace(0.0, S_LO, 8, endpoint=False),
                         np.linspace(S_LO, S_FULL, 14)])
camp = s_vals >= S_LO
HOST_BETA, HOST_5100 = 3.0, 1.825
HOST = HOST_5100 * (bands / 5100.0) ** HOST_BETA

WL = np.logspace(np.log10(800), np.log10(20000), 80)
nuWL = C / (WL * ANGSTROM)
nub = C / (bands * ANGSTROM)
i_hi = len(s_vals) - 1
i_lo = int(np.argmax(camp))

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

tls.config['disk']['nphi'] = NPHI
OBS = load_pg0844()


def make_disc(sigma_phi):
    """Pillared disc at one azimuthal width, with the same eta_LP rescaling and
    viscous floor as test_fluxflux_embedded.py. sigma_phi <= 0 gives the bowl."""
    dk = tls.build_disk(hlamp=HLAMP)
    if sigma_phi > 0:
        tls.add_pillars(dk, N_PILLARS, H_P, SIGMA_R, sigma_phi, r_min=5.0)
    dk.d = DL_MPC * PC_TO_LD
    grid = DiscGrid(dk)
    dk.tx_base *= (ETA_LP / eta_lp(grid)) ** 0.25
    dk.tv_base *= 1.5
    grid.tv2 *= 1.5
    return dk, grid


def sweep(grid, t4_spot, d0, tag):
    """Band fluxes over the lamp sweep, plus the bright-state temperature map.
    States are processed one at a time so the refined grid stays in memory."""
    dk = grid.dk
    tx0 = dk.tx_base.copy()
    diffuse = diffuse_flux(d0, bands)
    F = np.zeros((len(s_vals), len(bands)))
    T_bright = None
    for j, s in enumerate(s_vals):
        dk.tx_base = tx0 * s
        T = (dk.get_temperature(grid.r2, grid.p2, compute_shadows=True) ** 4
             + t4_spot) ** 0.25
        F[j] = grid.sed(T, bands) + diffuse
        if j == i_hi:
            T_bright = T
        print(f"  [{tag}] state {j+1}/{len(s_vals)}", flush=True)
    dk.tx_base = tx0
    return F, T_bright


def recover_host(F_states):
    """Add the true host, build the PyROA driver, fit the campaign states and
    extrapolate to the non-variability anchor."""
    Fh = F_states + HOST
    X = driving_lc(Fh, camp)
    fit = np.array([np.polyfit(X[camp], Fh[camp, k], 1)
                    for k in range(len(bands))])
    Xnv = nonvar_anchor(fit)
    return fit[:, 0] * Xnv + fit[:, 1], Fh, Xnv


# ---------------------------------------------------------------------------
# Run every sigma_phi, then the bowl baseline
# ---------------------------------------------------------------------------
res = {}
for sp in SIGMA_PHI_SCAN:
    print(f"\n=== sigma_phi = {sp} ===", flush=True)
    dk, grid = make_disc(sp)
    t4, areas, tpk = spot_t4(grid, M_EMBEDDED)
    d0 = diffuse_amplitude(grid, F_DIFF, DNU)
    F, T_bright = sweep(grid, t4, d0, f'sp={sp}')
    # shadow covering fraction at the nominal state, for the record
    smask = dk._compute_shadow_mask(grid.r2, grid.p2, grid.h2)
    res[sp] = {'dk': dk, 'grid': grid, 't4': t4, 'areas': areas, 'tpk': tpk,
               'F': F, 'T_bright': T_bright, 'd0': d0,
               'f_shadow': float(1.0 - smask.mean())}
    print(f"  A_bump = {areas.mean():.3e} cm^2, R_eff = "
          f"{np.sqrt(areas.mean()/np.pi)/AU_CM:.1f} AU, T_spot = "
          f"{tpk.min():.0f}-{tpk.max():.0f} K, shadow fraction = "
          f"{res[sp]['f_shadow']:.3f}", flush=True)

print("\n=== bowl (no pillars) ===", flush=True)
dk_b, grid_b = make_disc(0.0)
d0_b = diffuse_amplitude(grid_b, F_DIFF, DNU)
F_bowl, _ = sweep(grid_b, 0.0, d0_b, 'bowl')

# ---------------------------------------------------------------------------
# One global flux normalization, anchored on the reference sigma_phi so the
# embedded mass keeps the same physical meaning as in test_fluxflux_embedded.py
# ---------------------------------------------------------------------------
if OBS is None:
    BRIGHT_V = 9.0
else:
    BRIGHT_V = float(np.exp(np.interp(np.log(bands[iV]), np.log(OBS['lam']),
                                      np.log(OBS['bright']))))
FLUX_SCALE = (BRIGHT_V - HOST[iV]) / res[SIGMA_PHI_REF]['F'][:, iV].max()
M_PHYS = M_EMBEDDED * FLUX_SCALE
print(f"\nFLUX_SCALE = {FLUX_SCALE:.4g} (anchored at sigma_phi="
      f"{SIGMA_PHI_REF}); physical embedded mass = {M_PHYS:.3g} Msun per "
      f"pillar in total, i.e. {M_PHYS/M_STAR_UNIT:.0f} stars of "
      f"{M_STAR_UNIT:.0f} Msun")

for sp in SIGMA_PHI_SCAN:
    res[sp]['F'] *= FLUX_SCALE
    res[sp]['rec'], res[sp]['Fh'], res[sp]['Xnv'] = recover_host(res[sp]['F'])
F_bowl *= FLUX_SCALE
rec_bowl, Fh_bowl, Xnv_bowl = recover_host(F_bowl)

D_CM = res[SIGMA_PHI_REF]['dk'].d * LD_CM
NUL = 1e-26 * 4.0 * np.pi * D_CM**2

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
print(f"\nEmbedded mass held fixed at {M_PHYS:.3g} Msun per pillar "
      f"(L_pillar = {M_EMBEDDED*L_EDD_PER_MSUN*FLUX_SCALE:.2e} erg/s).")
print(f"{'sigma_phi':>10}{'arc@15ld':>10}{'A_bump':>12}{'R_eff[AU]':>11}"
      f"{'T_spot[K]':>16}{'f_shadow':>10}{'host@V':>9}{'host@i':>9}")
ki, kV = int(np.argmin(np.abs(bands - 7545))), iV
for sp in SIGMA_PHI_SCAN:
    r = res[sp]
    print(f"{sp:>10.2f}{15.0*sp:>10.2f}{r['areas'].mean():>12.3e}"
          f"{np.sqrt(r['areas'].mean()/np.pi)/AU_CM:>11.1f}"
          f"{r['tpk'].min():>8.0f}-{r['tpk'].max():<7.0f}"
          f"{r['f_shadow']:>10.3f}{r['rec'][kV]:>9.3f}{r['rec'][ki]:>9.3f}")
print(f"{'bowl':>10}{'-':>10}{'-':>12}{'-':>11}{'-':>16}{0.0:>10.3f}"
      f"{rec_bowl[kV]:>9.3f}{rec_bowl[ki]:>9.3f}")
if OBS is not None:
    hV = float(np.exp(np.interp(np.log(bands[kV]), np.log(OBS['lam']),
                                np.log(OBS['host']))))
    hi = float(np.exp(np.interp(np.log(bands[ki]), np.log(OBS['lam']),
                                np.log(OBS['host']))))
    print(f"{'PG 0844':>10}{'-':>10}{'-':>12}{'-':>11}{'-':>16}{'-':>10}"
          f"{hV:>9.3f}{hi:>9.3f}")

# ---------------------------------------------------------------------------
# Figure: recovered host and input spot spectrum for every sigma_phi
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 7))
cols = list(plt.cm.plasma(np.linspace(0.08, 0.78, len(SIGMA_PHI_SCAN))))
handles = []
ymark = []
for c, sp in zip(cols, SIGMA_PHI_SCAN):
    r = res[sp]
    spot_wl = r['grid'].sed(r['t4'] ** 0.25, WL) * FLUX_SCALE
    ok = r['rec'] > 0
    ax.plot(WL, nuWL*spot_wl*NUL, '-', color=c, lw=2.5, alpha=0.9)
    ax.plot(bands[ok], nub[ok]*r['rec'][ok]*NUL, '*', color=c, ms=16, mec='k',
            mew=0.8, ls='none', zorder=10)
    handles.append(Line2D([], [], color=c, marker='*', ms=14, mec='k', mew=0.8,
                          ls='-', lw=2.5,
                          label=rf'$\sigma_\phi = {sp:.2f}$'))
    ymark.append(nub[ok]*r['rec'][ok]*NUL)
ok_b = rec_bowl > 0
ax.plot(bands[ok_b], nub[ok_b]*rec_bowl[ok_b]*NUL, 'D', mfc='white',
        mec='dimgray', mew=2.0, ms=9, ls='none')
overlay_pg0844(ax, OBS, NUL, host_only=True)
if OBS is not None:
    ymark.append(OBS['nu'] * OBS['host'] * NUL)
keys = [Line2D([], [], color='k', ls='-', lw=2.5,
               label=r'$\rm spot~emission~(input)$'),
        Line2D([], [], color='k', marker='*', ms=14, mec='k', ls='none',
               label=r'$\rm host~(flux{-}flux~fit)$'),
        Line2D([], [], color='dimgray', marker='D', mfc='white', mew=2.0, ms=9,
               ls='none', label=r'$\rm host~(bowl~fit)$'),
        PG_KEYS[0]]
ax.set_xscale('log'); ax.set_yscale('log')
ymark = np.concatenate(ymark)
ax.set_ylim(3e-3 * ymark.max(), 8.0 * ymark.max())
ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=20)
ax.set_ylabel(r'$\nu L_\nu~\rm [erg~s^{-1}]$', fontsize=20)
leg = ax.legend(handles=handles, fontsize=12, frameon=True, loc='lower right',
                title=r'$\rm pillar~azimuthal~width$', title_fontsize=11)
ax.add_artist(leg)
ax.legend(handles=keys, fontsize=10, frameon=False, loc='upper left')
ticks(ax)
plt.tight_layout()
plt.savefig('plots/test_fluxflux_sigmaphi_host.png', dpi=200,
            bbox_inches='tight')
plt.close()
print("\nSaved plots/test_fluxflux_sigmaphi_host.png")

# ---------------------------------------------------------------------------
# Line-of-sight temperature maps, one per sigma_phi, on a common color scale
# ---------------------------------------------------------------------------
for sp in SIGMA_PHI_SCAN:
    r = res[sp]
    los_map(r['grid'], r['T_bright'],
            f'plots/test_fluxflux_sigmaphi_tmap_s{sp:.2f}.png',
            notes=[rf'${N_PILLARS}~{{\rm pillars}},~h_p = {H_P}~{{\rm ld}},~'
                   rf'\sigma_r = {SIGMA_R}~{{\rm ld}},~'
                   rf'\sigma_\phi = {sp:.2f}~({15.0*sp:.2f}~{{\rm ld~arc~at}}~'
                   rf'r = 15~{{\rm ld}})$',
                   r'$\rm embedded~stars:$ ' + mass_label(M_PHYS)
                   + r' $\rm per~pillar,$ '
                   + rf'$T_{{\rm spot}}^{{\rm peak}} = '
                     rf'{r["tpk"].min()/1e3:.1f}-{r["tpk"].max()/1e3:.1f}'
                     rf'~{{\rm kK}}$'],
            lim=20.0, vmin=5.0e2, vmax=4.0e4)
