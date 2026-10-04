"""Embedded-star flux-flux test with the spot size set by the luminosity.

test_fluxflux_embedded.py forces all of an embedded star's luminosity through
its own pillar's Gaussian footprint, which makes the spot temperature a free
knob: a narrower pillar gives a hotter spot at the same luminosity. That is not
what a radiation-hydrodynamic calculation would do. There the star is buried at
the midplane, its radiation diffuses out through the disc, and the size of the
warmed patch follows from the luminosity and the burial depth, not from the
shape of the bump.

This script replaces the spot prescription with

    sigma T_spot^4 = L h / (4 pi R^3),   R = sqrt(d^2 + h^2),

a point source at (r_p, phi_p, z=0) seen by a surface element at in-plane
separation d and height h. That integrates to exactly L/2 over the upper
surface, so energy conservation is built in rather than imposed, and the patch
width is no longer a free parameter. The consequence is that T_spot is pinned
near the local disc temperature: a brighter star warms a bigger patch out to the
same ambient floor rather than making a hotter spot.

Each mass is separately normalized so that its bright-state total matches the
observed PG 0844+349 bright state, so every run is a legitimate model of the
object and the comparison isolates the SHAPE of the recovered constant
component, which is the thing that is no longer free.

Outputs, in plots/ :
  test_fluxflux_rhd_host.png       recovered "host" for every mass and pillar
                                   width, against the PG 0844 host
  test_fluxflux_rhd_tmap_s<w>.png  line-of-sight temperature map
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pillardisk.pillar_disk import C, ANGSTROM, PC_TO_LD
import test_lag_spectrum as tls
from fluxflux_lib import (LD_CM, AU_CM, SIGMA_SB, L_EDD_PER_MSUN, DiscGrid,
                          spot_t4_rhd, heating_radius, diffuse_amplitude,
                          diffuse_flux, eta_lp, driving_lc, nonvar_anchor,
                          load_pg0844, overlay_pg0844, PG_KEYS, mass_label,
                          ticks, los_map)

# --- geometry ---------------------------------------------------------------
H_P, SIGMA_R = 0.5, 0.25
SIGMA_PHI_SCAN = [0.30, 0.05]       # wide (parent script) and RHD-preferred
HLAMP = 50 * 0.00399
N_PILLARS = 100
DL_MPC = 287.0
ETA_LP = 50.0
NPHI = 720
R_BELT = 15.0                        # light days, for the reported diagnostics

# --- embedded stars: TOTAL mass per pillar, in model (pre-scaling) units ----
M_SCAN = [5.0e3, 5.0e4, 2.0e5, 1.0e6]
M_FID = 2.0e5                        # used for the temperature maps
M_STAR_UNIT = 200.0

if '--quick' in sys.argv:
    SIGMA_PHI_SCAN, NPHI, M_SCAN = [0.30], 360, [5.0e4, 2.0e5]
    M_FID = 2.0e5
    print("QUICK MODE: coarse grid, one width, results are not final")

# --- bands, lamp sweep, diffuse glow, true host: as in the parent script ----
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

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

tls.config['disk']['nphi'] = NPHI
OBS = load_pg0844()
if OBS is None:
    BRIGHT_V = 9.0
else:
    BRIGHT_V = float(np.exp(np.interp(np.log(bands[iV]), np.log(OBS['lam']),
                                      np.log(OBS['bright']))))


def make_disc(sigma_phi):
    dk = tls.build_disk(hlamp=HLAMP)
    if sigma_phi > 0:
        tls.add_pillars(dk, N_PILLARS, H_P, SIGMA_R, sigma_phi, r_min=5.0)
    dk.d = DL_MPC * PC_TO_LD
    grid = DiscGrid(dk)
    dk.tx_base *= (ETA_LP / eta_lp(grid)) ** 0.25
    dk.tv_base *= 1.5
    grid.tv2 *= 1.5
    return dk, grid


def lamp_T4(grid, tag):
    """T^4 of the viscous + irradiated disc for every lamp state, kept in full
    because the spot term is added afterwards for each mass."""
    dk = grid.dk
    tx0 = dk.tx_base.copy()
    out = np.empty((len(s_vals),) + grid.r2.shape)
    for j, s in enumerate(s_vals):
        dk.tx_base = tx0 * s
        out[j] = dk.get_temperature(grid.r2, grid.p2,
                                    compute_shadows=True) ** 4
        print(f"  [{tag}] state {j+1}/{len(s_vals)}", flush=True)
    dk.tx_base = tx0
    return out


def recover_host(F_states):
    Fh = F_states + HOST
    X = driving_lc(Fh, camp)
    fit = np.array([np.polyfit(X[camp], Fh[camp, k], 1)
                    for k in range(len(bands))])
    Xnv = nonvar_anchor(fit)
    return fit[:, 0] * Xnv + fit[:, 1], Fh


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
res = {}
for sp in SIGMA_PHI_SCAN:
    print(f"\n=== sigma_phi = {sp} ===", flush=True)
    dk, grid = make_disc(sp)
    T4 = lamp_T4(grid, f'sp={sp}')
    d0 = diffuse_amplitude(grid, F_DIFF, DNU)
    diffuse = diffuse_flux(d0, bands)
    # the spot profile shape does not depend on the mass, so build it once at a
    # reference mass and scale (F is linear in L)
    t4_ref, caught, tpk_ref = spot_t4_rhd(grid, M_SCAN[0])
    ir = int(np.argmin(np.abs(dk.r - R_BELT)))
    T_amb = float(np.mean(T4[i_hi][ir] ** 0.25))   # azimuthal mean at the belt
    print(f"  discrete profile captured {caught.mean()*100:.1f}% of L/2 before "
          f"renormalization (min {caught.min()*100:.1f}%, "
          f"max {caught.max()*100:.1f}%)")
    print(f"  ambient T at r={R_BELT} ld: {T_amb:.0f} K", flush=True)
    res[sp] = {'dk': dk, 'grid': grid, 'T4': T4, 'diffuse': diffuse,
               't4_ref': t4_ref, 'T_amb': T_amb, 'runs': {}}

    for m in M_SCAN:
        t4_spot = t4_ref * (m / M_SCAN[0])
        F = np.array([grid.sed((T4[j] + t4_spot) ** 0.25, bands) + diffuse
                      for j in range(len(s_vals))])
        # normalize this run on its own so its bright state matches PG 0844
        fs = (BRIGHT_V - HOST[iV]) / F[:, iV].max()
        rec, Fh = recover_host(F * fs)
        t_pk = float((t4_spot.max()) ** 0.25)
        r_heat = heating_radius(m * fs, T_amb)
        res[sp]['runs'][m] = {'F': F * fs, 'rec': rec, 'Fh': Fh, 'fs': fs,
                              'm_phys': m * fs, 't_pk': t_pk,
                              'r_heat': r_heat, 't4_spot': t4_spot}
        print(f"  M_model={m:.3g} -> FLUX_SCALE={fs:.4g}, M_phys="
              f"{m*fs:.3g} Msun/pillar, T_spot_peak={t_pk:.0f} K, "
              f"R_heat={r_heat/AU_CM:.0f} AU = {r_heat/LD_CM:.2f} ld",
              flush=True)

D_CM = res[SIGMA_PHI_SCAN[0]]['dk'].d * LD_CM
NUL = 1e-26 * 4.0 * np.pi * D_CM**2

# ---------------------------------------------------------------------------
# Report: recovered host against the observed one
# ---------------------------------------------------------------------------
kV, ki = iV, int(np.argmin(np.abs(bands - 7545)))
if OBS is not None:
    hV = float(np.exp(np.interp(np.log(bands[kV]), np.log(OBS['lam']),
                                np.log(OBS['host']))))
    hi = float(np.exp(np.interp(np.log(bands[ki]), np.log(OBS['lam']),
                                np.log(OBS['host']))))
    hU = float(np.exp(np.interp(np.log(bands[3]), np.log(OBS['lam']),
                                np.log(OBS['host']))))
print("\nRecovered constant ('host'), mJy. Each run is separately normalized so "
      "its\nbright state matches PG 0844, so this compares the SHAPE.")
print(f"{'sigma_phi':>10}{'M_phys':>11}{'T_peak[K]':>11}{'R_heat[ld]':>12}"
      f"{'host@U':>9}{'host@V':>9}{'host@i':>9}{'U/i':>7}")
for sp in SIGMA_PHI_SCAN:
    for m in M_SCAN:
        r = res[sp]['runs'][m]
        print(f"{sp:>10.2f}{r['m_phys']:>11.3g}{r['t_pk']:>11.0f}"
              f"{r['r_heat']/LD_CM:>12.2f}{r['rec'][3]:>9.3f}"
              f"{r['rec'][kV]:>9.3f}{r['rec'][ki]:>9.3f}"
              f"{r['rec'][3]/max(r['rec'][ki], 1e-9):>7.3f}")
if OBS is not None:
    print(f"{'PG 0844':>10}{'-':>11}{'-':>11}{'-':>12}{hU:>9.3f}{hV:>9.3f}"
          f"{hi:>9.3f}{hU/hi:>7.3f}")

# ---------------------------------------------------------------------------
# Figure: recovered host and input spot spectrum
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 7))
cols = list(plt.cm.viridis(np.linspace(0.1, 0.85, len(M_SCAN))))
lstyle = {SIGMA_PHI_SCAN[0]: '-'}
for k, sp in enumerate(SIGMA_PHI_SCAN[1:], start=1):
    lstyle[sp] = ['--', ':', '-.'][k - 1]
handles, ymark = [], []
for sp in SIGMA_PHI_SCAN:
    for c, m in zip(cols, M_SCAN):
        r = res[sp]['runs'][m]
        spot_wl = res[sp]['grid'].sed(r['t4_spot'] ** 0.25, WL) * r['fs']
        ok = r['rec'] > 0
        ax.plot(WL, nuWL*spot_wl*NUL, lstyle[sp], color=c, lw=2.2, alpha=0.85)
        ax.plot(bands[ok], nub[ok]*r['rec'][ok]*NUL,
                '*' if sp == SIGMA_PHI_SCAN[0] else 'P', color=c, ms=15,
                mec='k', mew=0.8, ls='none', zorder=10)
        ymark.append(nub[ok]*r['rec'][ok]*NUL)
for c, m in zip(cols, M_SCAN):
    handles.append(Line2D([], [], color=c, lw=2.5, marker='*', ms=14, mec='k',
                          mew=0.8,
                          label=mass_label(res[SIGMA_PHI_SCAN[0]]['runs'][m]['m_phys'])))
overlay_pg0844(ax, OBS, NUL, host_only=True)
if OBS is not None:
    ymark.append(OBS['nu'] * OBS['host'] * NUL)
keys = [Line2D([], [], color='k', ls='-', lw=2.2,
               label=r'$\rm spot~emission~(input)$'),
        Line2D([], [], color='k', marker='*', ms=14, mec='k', ls='none',
               label=rf'$\rm host~fit,~\sigma_\phi={SIGMA_PHI_SCAN[0]:.2f}$'),
        PG_KEYS[0]]
if len(SIGMA_PHI_SCAN) > 1:
    keys.insert(2, Line2D([], [], color='k', marker='P', ms=13, mec='k',
                          ls='none',
                          label=rf'$\rm host~fit,~\sigma_\phi='
                                rf'{SIGMA_PHI_SCAN[1]:.2f}$'))
ax.set_xscale('log'); ax.set_yscale('log')
ymark = np.concatenate(ymark)
ax.set_ylim(3e-3 * ymark.max(), 8.0 * ymark.max())
ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=20)
ax.set_ylabel(r'$\nu L_\nu~\rm [erg~s^{-1}]$', fontsize=20)
leg = ax.legend(handles=handles, fontsize=11, frameon=True, loc='lower right',
                title=r'$\rm embedded~object~mass~per~pillar$',
                title_fontsize=10)
ax.add_artist(leg)
ax.legend(handles=keys, fontsize=10, frameon=False, loc='upper left')
ticks(ax)
plt.tight_layout()
plt.savefig('plots/test_fluxflux_rhd_host.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_fluxflux_rhd_host.png")

# ---------------------------------------------------------------------------
# Temperature maps at the fiducial mass
# ---------------------------------------------------------------------------
for sp in SIGMA_PHI_SCAN:
    r = res[sp]['runs'][M_FID]
    T_map = (res[sp]['T4'][i_hi] + r['t4_spot']) ** 0.25
    los_map(res[sp]['grid'], T_map,
            f'plots/test_fluxflux_rhd_tmap_s{sp:.2f}.png',
            notes=[rf'${N_PILLARS}~{{\rm pillars}},~h_p = {H_P}~{{\rm ld}},~'
                   rf'\sigma_r = {SIGMA_R}~{{\rm ld}},~'
                   rf'\sigma_\phi = {sp:.2f}$',
                   r'$\rm buried~sources:$ ' + mass_label(r['m_phys'])
                   + r' $\rm per~pillar,$ '
                   + rf'$R_{{\rm heat}} = {r["r_heat"]/LD_CM:.2f}~{{\rm ld}},~'
                     rf'T_{{\rm spot}}^{{\rm peak}} = {r["t_pk"]/1e3:.1f}'
                     rf'~{{\rm kK}}$'],
            lim=20.0, vmin=5.0e2, vmax=4.0e4)
