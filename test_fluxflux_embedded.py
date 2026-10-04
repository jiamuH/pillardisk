"""Flux-flux decomposition when the pillars carry embedded-star hot spots.

Physical picture: stars embedded in the AGN disc accrete, and their luminosity
has to go somewhere even when the star is fully buried in dense gas. Here it is
thermalized and re-radiated by the pillar bump that sits above the star. Each
pillar therefore carries a constant, lamp-independent blackbody component on top
of the viscous plus irradiated disc emission,

    sigma * T_spot^4 (r, phi) = sum_p  L_p * g_p(r, phi) / A_p ,
    A_p = Integral of g_p over the disc surface (the pillar's own Gaussian
          footprint area, in cm^2),

with L_p = M_embedded * 1.26e38 erg/s, the Eddington luminosity of the embedded
stellar mass assigned to that pillar. The pillar surface radiates the sum of all
heating channels, so temperatures add as
    T^4 = T_visc^4 + T_irr^4 + T_spot^4 ,
which also means the spot suppresses the fractional response of the pillar
surface to the varying lamp.

Because the spot component never varies with the lamp, a PyROA-style flux-flux
extrapolation dumps it entirely into the recovered constant ("host") SED. The
recovered host then comes out blackbody-like and far brighter than the true
host, which is the behavior seen in the PG 0844+349 flux-flux analysis.

Everything else (disc, pillars, lamp sweep, driving light curve, non-variability
anchor, true host SED) is identical to test_fluxflux_mock.py, so the figures can
be compared one to one with plots/test_fluxflux_sed.png.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pillardisk.pillar_disk import C, ANGSTROM, PC_TO_LD
from test_lag_spectrum import build_disk, add_pillars
from fluxflux_lib import (LD_CM, SIGMA_SB, L_EDD_PER_MSUN, DiscGrid, spot_t4,
                          diffuse_amplitude, diffuse_flux, lamp_states, eta_lp,
                          driving_lc as _driving_lc, nonvar_anchor, load_pg0844,
                          overlay_pg0844 as _overlay_pg0844, PG_KEYS,
                          mass_label, ticks as _ticks, los_map as _los_map)

# --- disc / pillar setup, identical to test_fluxflux_mock.py ---
H_P, SIGMA_R, SIGMA_PHI = 0.5, 0.25, 0.3
HLAMP = 50 * 0.00399           # lamp height = 50 r_g (r_g = 0.00399 ld, 7e7 Msun)
N_PILLARS = 100
Z_PG0844, DL_MPC = 0.064, 287.0
ETA_LP = 50.0                  # bolometric lamp / viscous power at the nominal state

# --- embedded stars -------------------------------------------------------
# Mass of the SINGLE embedded object under each pillar, radiating at its
# Eddington luminosity and thermalized by the bump. One pillar is raised by one
# object, so this is an individual mass. Only the product
# L_pillar = 1.26e38 (M/Msun) erg/s enters the model. Note that the values the
# flux-flux match demands are far above any stellar mass: a real star of at most
# a few hundred Msun supplies of order 1% of the observed constant component
# (see debug_one_star_per_pillar.py), so the upper end of this scan corresponds
# to intermediate-mass black holes, not stars. The scan starts from zero, the
# shadow-only case of test_fluxflux_mock.py.
M_STAR_UNIT = 200.0            # a massive star, for comparison against the scan
M_EMBEDDED_SCAN = [0.0, 2.0e2, 5.0e3, 5.0e4, 2.0e5]   # Msun per pillar
M_EMBEDDED_FID = 5.0e4                                 # Msun per pillar, fiducial

# --- observed bands (PG 0844-like, rest-frame effective wavelengths, A) ---
bands = np.array([1928.0, 2246.0, 2600.0, 3465.0, 4392.0,
                  5468.0, 6215.0, 7545.0, 8700.0])
band_lab = [r'$\rm UVW2$', r'$\rm UVM2$', r'$\rm UVW1$', r'$U$', r'$B$',
            r'$V$', r'$r$', r'$i$', r'$z$']
band_col = list(plt.cm.turbo(np.linspace(0.06, 0.95, len(bands))))
iV = int(np.argmin(np.abs(bands - 5100.0)))

# --- diffuse shadow glow, identical to test_fluxflux_mock.py ---
F_DIFF = 0.3
DNU = C / (bands.min() * ANGSTROM) - C / (bands.max() * ANGSTROM)

# --- lamp sweep and campaign window, identical to test_fluxflux_mock.py ---
S_FULL = 1.35
CAMP_FRAC = 0.20
S_HI = S_FULL
S_LO = S_FULL * (1.0 - CAMP_FRAC)
s_vals = np.concatenate([np.linspace(0.0, S_LO, 8, endpoint=False),
                         np.linspace(S_LO, S_HI, 14)])
CAMP_MIN = S_LO
camp = s_vals >= CAMP_MIN

# --- true host SED, identical to test_fluxflux_mock.py ---
HOST_BETA, HOST_5100 = 3.0, 1.825
HOST = HOST_5100 * (bands / 5100.0) ** HOST_BETA


def driving_lc(F):
    """PyROA driver over the campaign window defined in this script."""
    return _driving_lc(F, camp)


def overlay_pg0844(ax, o):
    """PG 0844 overlay at this script's nu L_nu conversion."""
    return _overlay_pg0844(ax, o, NUL)


# ---------------------------------------------------------------------------
# Build the pillared disc and the bowl baseline
# ---------------------------------------------------------------------------
disk = build_disk(hlamp=HLAMP)
add_pillars(disk, N_PILLARS, H_P, SIGMA_R, SIGMA_PHI, r_min=5.0)
disk.d = DL_MPC * PC_TO_LD
bowl = build_disk(hlamp=HLAMP)
bowl.d = DL_MPC * PC_TO_LD
print(f"Disc: rin={disk.rin:.3f} rout={disk.rout} ld, nr={disk.nr} "
      f"nphi={disk.nphi}, {len(disk.pillars)} pillars; D_L={DL_MPC} Mpc")

gp = DiscGrid(disk)
gb = DiscGrid(bowl)

eta0 = eta_lp(gp)
f_eta = (ETA_LP / eta0) ** 0.25
disk.tx_base *= f_eta
bowl.tx_base *= f_eta
disk.tv_base *= 1.5
bowl.tv_base *= 1.5
gp.tv2 *= 1.5
gb.tv2 *= 1.5
print(f"eta_LP: original {eta0:.1f} -> rescaled by {f_eta:.3f} on T_irr; "
      f"viscous floor raised x1.5 -> lamp/viscous = {eta_lp(gp):.1f}")

# ---------------------------------------------------------------------------
# Lamp states (embedded-star independent) and the constant diffuse glow
# ---------------------------------------------------------------------------
print("Computing lamp states for the pillared disc...", flush=True)
T4_p = lamp_states(gp, s_vals, 'pillars')
print("Computing lamp states for the bowl...", flush=True)
T4_b = lamp_states(gb, s_vals, 'bowl')
D0_p = diffuse_amplitude(gp, F_DIFF, DNU)
D0_b = diffuse_amplitude(gb, F_DIFF, DNU)
diff_p = diffuse_flux(D0_p, bands)
diff_b = diffuse_flux(D0_b, bands)

# ---------------------------------------------------------------------------
# Embedded-star spot maps
# ---------------------------------------------------------------------------
spot_maps = {}
print("\nEmbedded-star hot spots (Eddington luminosity, thermalized by the bump).")
print("M_emb is the mass of the single embedded object under each pillar; "
      f"N_star is\nhow many {M_STAR_UNIT:.0f} Msun stars would be needed to "
      "match that luminosity,\nwhich is only a reference, not the model setup.")
print(f"{'M_emb/pillar':>14}{'N_star':>9}{'L_pillar':>12}{'L_total':>12}"
      f"{'A_bump [cm^2]':>16}{'R_eff [AU]':>12}{'T_spot [K]':>14}")
for m_emb in M_EMBEDDED_SCAN:
    t4, areas, tpk = spot_t4(gp, m_emb)
    spot_maps[m_emb] = t4
    r_eff_au = np.sqrt(areas.mean() / np.pi) / 1.496e13
    print(f"{m_emb:>14.3g}{m_emb/M_STAR_UNIT:>9.0f}{m_emb*L_EDD_PER_MSUN:>12.2e}"
          f"{m_emb*L_EDD_PER_MSUN*N_PILLARS:>12.2e}{areas.mean():>16.3e}"
          f"{r_eff_au:>12.1f}{tpk.min():>7.0f}-{tpk.max():<6.0f}")

# ---------------------------------------------------------------------------
# Flux-flux sweep for every embedded mass (pillars) and for the bowl (no spots)
# ---------------------------------------------------------------------------
def band_fluxes(grid, T4_lamp, t4_spot, diffuse):
    """Band fluxes (n_state, n_band) plus the lamp-off constant SED."""
    F = np.array([grid.sed((T4 + t4_spot) ** 0.25, bands) + diffuse
                  for T4 in T4_lamp])
    F0 = grid.sed((grid.tv2**4 + t4_spot) ** 0.25, bands) + diffuse
    return F, F0


F_bowl, F0_bowl = band_fluxes(gb, T4_b, 0.0, diff_b)
runs = {}
for m_emb in M_EMBEDDED_SCAN:
    runs[m_emb] = band_fluxes(gp, T4_p, spot_maps[m_emb], diff_p)
    print(f"flux-flux sweep done for M_embedded = {m_emb:.3g} Msun/pillar",
          flush=True)

# Absolute scale: the model SED normalization is arbitrary, so rescale so the
# bright-state model total (disc + spots + true host) matches PG 0844's observed
# bright-state flux at the model's V wavelength. FLUX_SCALE is one global
# constant applied to every run, so the underlying disc stays identical across
# the embedded mass scan and only the spot term differs between runs. The same
# factor multiplies every luminosity in the model, so the physical embedded mass
# implied by a run is M_embedded * FLUX_SCALE.
OBS = load_pg0844()
if OBS is None:
    BRIGHT_V = 9.0                      # PG 0844 bright-state V, mJy (fallback)
else:
    BRIGHT_V = float(np.exp(np.interp(np.log(bands[iV]), np.log(OBS['lam']),
                                      np.log(OBS['bright']))))
FLUX_SCALE = (BRIGHT_V - HOST[iV]) / runs[M_EMBEDDED_FID][0][:, iV].max()
print(f"\nAnchor: observed bright-state flux at {bands[iV]:.0f} A is "
      f"{BRIGHT_V:.2f} mJy; model true host there is {HOST[iV]:.2f} mJy")
print(f"\nFLUX_SCALE = {FLUX_SCALE:.3g}  (implied physical embedded mass = "
      f"model mass x {FLUX_SCALE:.3g})")
F_bowl *= FLUX_SCALE
F0_bowl *= FLUX_SCALE
runs = {m: (F * FLUX_SCALE, F0 * FLUX_SCALE) for m, (F, F0) in runs.items()}


def recover_host(F_states):
    """Add the true host, build the PyROA driver, fit the campaign states and
    extrapolate to the non-variability anchor. Returns the recovered constant
    SED, the driver, the fits and the anchor."""
    Fh = F_states + HOST
    X = driving_lc(Fh)
    fit = np.array([np.polyfit(X[camp], Fh[camp, k], 1)
                    for k in range(len(bands))])
    Xnv = nonvar_anchor(fit)
    return fit[:, 0] * Xnv + fit[:, 1], X, fit, Xnv, Fh


rec_bowl, X_bowl, fit_bowl, Xnv_bowl, Fh_bowl = recover_host(F_bowl)
rec = {}
for m_emb in M_EMBEDDED_SCAN:
    rec[m_emb] = recover_host(runs[m_emb][0])

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
rec_nospot = rec[0.0][0]
print("\nRecovered constant ('host') SED vs the true host, mJy. 'excess' is the "
      "extra recovered flux over the spots-off (shadow-only) case:")
for m_emb in M_EMBEDDED_SCAN:
    r_h = rec[m_emb][0]
    print(f"  M_embedded = {m_emb:.3g} Msun/pillar "
          f"(physical {m_emb*FLUX_SCALE:.3g}), X_non-var = {rec[m_emb][3]:.2f}")
    for k, lam in enumerate(bands):
        print(f"    {lam:6.0f} A : recovered={r_h[k]:7.3f}  "
              f"true_host={HOST[k]:6.3f}  ratio={r_h[k]/HOST[k]:6.2f}  "
              f"excess={r_h[k]-rec_nospot[k]:7.3f}")
print("\nBowl (no pillars, no spots) recovered host, mJy:")
for k, lam in enumerate(bands):
    print(f"  {lam:6.0f} A : recovered={rec_bowl[k]:7.3f}  "
          f"true_host={HOST[k]:6.3f}")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

WL = np.logspace(np.log10(800), np.log10(20000), 80)
nuWL = C / (WL * ANGSTROM)
WL_host = np.logspace(np.log10(1000), np.log10(20000), 70)
nuWL_host = C / (WL_host * ANGSTROM)
nub = C / (bands * ANGSTROM)
D_CM = disk.d * LD_CM
NUL = 1e-26 * 4.0 * np.pi * D_CM**2          # mJy -> erg s^-1 with the nu factor
HOST_WL = HOST_5100 * (WL_host / 5100.0) ** HOST_BETA

i_hi = len(s_vals) - 1                        # brightest state
i_lo = int(np.argmax(camp))                   # faintest observed state

def overlay_pg0844(ax, o):
    """Overplot the observed PG 0844 decomposition in nu L_nu. The raw bright
    and faint states are drawn faded; the host-subtracted states and the
    flux-flux host are opaque."""
    if o is None:
        return
    nl = o['nu'] * NUL
    for key, col in [('bright', 'royalblue'), ('faint', 'navy')]:
        ax.errorbar(o['lam'], nl*o[key], yerr=nl*o[key + '_e'], fmt='o',
                    mfc='white', mec=col, ms=7, mew=1.8, ecolor=col,
                    alpha=0.3, ls='none', zorder=8)
        sub = o[key] - o['host']
        err = np.hypot(o[key + '_e'], o['host_e'])
        ax.errorbar(o['lam'], nl*sub, yerr=nl*err, fmt='o', mfc='white',
                    mec=col, ms=7, mew=1.8, ecolor=col, ls='none', zorder=9)
    ax.errorbar(o['lam'], nl*o['host'], yerr=nl*o['host_e'], fmt='s',
                color='darkred', mec='k', mew=0.8, ms=8, ecolor='darkred',
                ls='none', zorder=11)


def sed_curve(grid, T4_lamp_state, t4_spot, d0, wl):
    """Full model SED (mJy) on a wavelength grid, at the plot normalization."""
    return (grid.sed((T4_lamp_state + t4_spot) ** 0.25, wl)
            + diffuse_flux(d0, wl)) * FLUX_SCALE


def spot_only_curve(grid, t4_spot, wl):
    """The embedded-star spot emission alone (mJy), a superposition of the
    per-pillar blackbodies with no viscous or irradiation heating."""
    return grid.sed(t4_spot ** 0.25, wl) * FLUX_SCALE


# --- Figure 1: SED view at the fiducial embedded mass ---------------------
t4f = spot_maps[M_EMBEDDED_FID]
rec_p, X_p, fit_p, Xnv_p, Fh_p = rec[M_EMBEDDED_FID]
F_p = runs[M_EMBEDDED_FID][0]

TpB = sed_curve(gp, T4_p[i_hi], t4f, D0_p, WL)
TpF = sed_curve(gp, T4_p[i_lo], t4f, D0_p, WL)
TbB = sed_curve(gb, T4_b[i_hi], 0.0, D0_b, WL)
TbF = sed_curve(gb, T4_b[i_lo], 0.0, D0_b, WL)
SPOT_WL = spot_only_curve(gp, t4f, WL)
inf_pB = Fh_p[i_hi] - rec_p
inf_pF = Fh_p[i_lo] - rec_p
inf_bB = Fh_bowl[i_hi] - rec_bowl
inf_bF = Fh_bowl[i_lo] - rec_bowl

fig, axs = plt.subplots(figsize=(8.5, 7))
axs.plot(WL, nuWL*TpB*NUL, '-', color='royalblue', lw=2.5,
         label=r'$\rm high~AGN~(pillars+spots)$')
axs.plot(WL, nuWL*TpF*NUL, '-', color='navy', lw=2.5,
         label=r'$\rm low~AGN~(pillars+spots)$')
axs.plot(WL, nuWL*TbB*NUL, '--', color='royalblue', lw=2.0, alpha=0.6,
         label=r'$\rm high~AGN~(bowl)$')
axs.plot(WL, nuWL*TbF*NUL, '--', color='navy', lw=2.0, alpha=0.6,
         label=r'$\rm low~AGN~(bowl)$')
axs.plot(WL, nuWL*SPOT_WL*NUL, '-', color='seagreen', lw=3.0, alpha=0.9,
         label=r'$\rm embedded{-}star~spots~(input)$')
axs.plot(WL_host, nuWL_host*HOST_WL*NUL, '-', color='firebrick', lw=2.5,
         label=r'$\rm host~(true)$')
mp, mb = rec_p > 0, rec_bowl > 0
axs.plot(bands[mp], nub[mp]*rec_p[mp]*NUL, '*', color='orange', ms=16, mec='k',
         mew=0.8, ls='none', zorder=12, label=r'$\rm host~(pillar+spot~fit)$')
axs.plot(bands[mb], nub[mb]*rec_bowl[mb]*NUL, 'D', mfc='white',
         mec='darkorange', mew=2.0, ms=9, ls='none',
         label=r'$\rm host~(bowl~fit)$')
axs.plot(bands, nub*inf_pB*NUL, 'x', color='royalblue', ms=9, mew=2.2, ls=':', lw=1.3)
axs.plot(bands, nub*inf_pF*NUL, 'x', color='navy', ms=9, mew=2.2, ls=':', lw=1.3)
axs.plot(bands, nub*inf_bB*NUL, '+', color='royalblue', ms=11, mew=2.2, ls=':', lw=1.3)
axs.plot(bands, nub*inf_bF*NUL, '+', color='navy', ms=11, mew=2.2, ls=':', lw=1.3)
overlay_pg0844(axs, OBS)
inf_keys = [Line2D([], [], color='k', marker='x', ms=9, mew=2.2, ls='none',
                   label=r'$\rm inferred~AGN~(pillars)$'),
            Line2D([], [], color='k', marker='+', ms=11, mew=2.2, ls='none',
                   label=r'$\rm inferred~AGN~(bowl)$')] + PG_KEYS
axs.set_xscale('log'); axs.set_yscale('log')
allnu = np.concatenate([nuWL*TpB, nuWL*TbF, nuWL_host*HOST_WL]) * NUL
axs.set_ylim(0.4 * allnu.min(), 3.0 * allnu.max())
axs.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=20)
axs.set_ylabel(r'$\nu L_\nu~\rm [erg~s^{-1}]$', fontsize=20)
leg1 = axs.legend(fontsize=10, frameon=True, loc='upper right', ncol=2)
axs.add_artist(leg1)
axs.legend(handles=inf_keys, fontsize=9, frameon=False, loc='lower left')
_ticks(axs)
plt.tight_layout()
plt.savefig('plots/test_fluxflux_embedded_sed.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_fluxflux_embedded_sed.png")

# --- Figure 2: recovered host vs embedded mass ----------------------------
fig2, ax2 = plt.subplots(figsize=(8.5, 7))
mass_col = list(plt.cm.viridis(np.linspace(0.1, 0.85, len(M_EMBEDDED_SCAN))))
ymark = []
mass_handles = []
for c, m_emb in zip(mass_col, M_EMBEDDED_SCAN):
    r_h = rec[m_emb][0]
    m_ok = r_h > 0
    if m_emb == 0.0:
        lab = r'$\rm none~(shadows~only)$'
    else:
        lab = mass_label(m_emb * FLUX_SCALE)
        sp_wl = spot_only_curve(gp, spot_maps[m_emb], WL)
        ax2.plot(WL, nuWL*sp_wl*NUL, '-', color=c, lw=2.5, alpha=0.9)
    ax2.plot(bands[m_ok], nub[m_ok]*r_h[m_ok]*NUL, '*', color=c, ms=16,
             mec='k', mew=0.8, ls='none', zorder=10)
    mass_handles.append(Line2D([], [], color=c, marker='*', ms=14, mec='k',
                               mew=0.8, ls='none', label=lab))
    ymark.append(nub[m_ok]*r_h[m_ok]*NUL)
ax2.plot(WL_host, nuWL_host*HOST_WL*NUL, '-', color='firebrick', lw=2.5)
ax2.plot(bands, nub*rec_bowl*NUL, 'D', mfc='white', mec='dimgray', mew=2.0,
         ms=9, ls='none')
if OBS is not None:
    overlay_pg0844(ax2, OBS)
    ymark.append(OBS['nu'] * OBS['bright'] * NUL)
key2 = [Line2D([], [], color='k', ls='-', lw=2.5,
               label=r'$\rm spot~emission~(input)$'),
        Line2D([], [], color='k', marker='*', ms=14, mec='k', ls='none',
               label=r'$\rm host~(flux{-}flux~fit)$'),
        Line2D([], [], color='firebrick', ls='-', lw=2.5,
               label=r'$\rm host~(true)$'),
        Line2D([], [], color='dimgray', marker='D', mfc='white', mew=2.0,
               ms=9, ls='none', label=r'$\rm host~(bowl~fit)$')] + PG_KEYS
ax2.set_xscale('log'); ax2.set_yscale('log')
# y range set by the true host at the blue end and the brightest recovered
# point at the top; the faint spot curves run far below and are simply clipped
ymark = np.concatenate(ymark)
ax2.set_ylim(0.2 * (nuWL_host * HOST_WL * NUL).min(), 8.0 * ymark.max())
ax2.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=20)
ax2.set_ylabel(r'$\nu L_\nu~\rm [erg~s^{-1}]$', fontsize=20)
leg2 = ax2.legend(handles=mass_handles, fontsize=11, frameon=True,
                  loc='lower right',
                  title=r'$\rm total~embedded~stellar~mass~per~pillar$'
                        '\n'
                        r'$\rm (summed~over~the~whole~population,~at~Eddington)$',
                  title_fontsize=9)
ax2.add_artist(leg2)
ax2.legend(handles=key2, fontsize=10, frameon=False, loc='upper left', ncol=2)
_ticks(ax2)
plt.tight_layout()
plt.savefig('plots/test_fluxflux_embedded_mass.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved plots/test_fluxflux_embedded_mass.png")

# --- Figure 3: the flux-flux diagram at the fiducial embedded mass ---------
Xf_p, Xb_p = X_p[camp].min(), X_p[camp].max()
Xf_b, Xb_b = X_bowl[camp].min(), X_bowl[camp].max()
fig3, axr = plt.subplots(figsize=(8, 7))
for k in range(len(bands)):
    c = band_col[k]
    axr.plot(X_p, Fh_p[:, k], 'o-', color=c, lw=2.5, ms=6, alpha=0.22)
    axr.plot(X_p[camp], Fh_p[camp, k], 'o-', color=c, lw=2.5, ms=6, alpha=0.9)
    axr.plot(X_bowl, Fh_bowl[:, k], marker='s', ls='--', color=c, mfc='white',
             lw=2.0, ms=5, alpha=0.18)
    axr.plot(X_bowl[camp], Fh_bowl[camp, k], marker='s', ls='--', color=c,
             mfc='white', lw=2.0, ms=5, alpha=0.6)
    xin = np.array([Xf_p, Xb_p])
    axr.plot(xin, fit_p[k, 0]*xin + fit_p[k, 1], ':', color=c, lw=1.5, alpha=0.7)
    axr.plot([Xnv_p, Xf_p], fit_p[k, 0]*np.array([Xnv_p, Xf_p]) + fit_p[k, 1],
             ':', color=c, lw=1.2, alpha=0.3)
    xib = np.array([Xf_b, Xb_b])
    axr.plot(xib, fit_bowl[k, 0]*xib + fit_bowl[k, 1], '-.', color=c, lw=1.3,
             alpha=0.6)
    axr.plot([Xnv_bowl, Xf_b],
             fit_bowl[k, 0]*np.array([Xnv_bowl, Xf_b]) + fit_bowl[k, 1],
             '-.', color=c, lw=1.0, alpha=0.25)
    axr.plot(Xnv_p, rec_p[k], marker='*', color=c, ms=16, mec='k', mew=0.8,
             ls='none', zorder=5)
    axr.plot(Xnv_bowl, rec_bowl[k], marker='o', mfc='white', mec=c, mew=2.0,
             ms=10, ls='none', zorder=5)
    axr.plot(Xnv_p, HOST[k], marker='D', color=c, ms=8, mec='k', mew=0.8,
             ls='none', zorder=6)
axr.axhline(0.0, color='gray', ls='--', lw=1.5, alpha=0.7)
axr.axvline(Xnv_p, color='gray', ls='--', lw=1.5, alpha=0.7)
axr.axvline(Xnv_bowl, color='gray', ls='-.', lw=1.2, alpha=0.6)
ytop = axr.get_ylim()[1]
for xv, lab in [(Xf_p, r'$X_{\rm faint}$'), (Xb_p, r'$X_{\rm bright}$')]:
    axr.axvline(xv, color='k', ls=':', lw=1.5, alpha=0.6)
    axr.text(xv, ytop, lab, rotation=90, va='top', ha='right', fontsize=11,
             color='k')
real_keys = [Line2D([], [], color='k', ls='-', marker='o',
                    label=r'$\rm pillars+spots$'),
             Line2D([], [], color='k', ls='--', marker='s', mfc='white',
                    label=r'$\rm bowl$'),
             Line2D([], [], color='k', marker='*', mec='k', ms=14, ls='none',
                    label=r'$\rm recovered,~pillar$'),
             Line2D([], [], color='k', marker='o', mfc='white', mec='k', ms=9,
                    ls='none', label=r'$\rm recovered,~bowl$'),
             Line2D([], [], color='k', marker='D', mec='k', ms=8, ls='none',
                    label=r'$\rm true~host$')]
band_handles = [Line2D([], [], color=band_col[k], marker='o', ls='-',
                       label=band_lab[k]) for k in range(len(bands))][::-1]
leg = axr.legend(handles=band_handles, fontsize=12, frameon=True, ncol=1,
                 loc='upper left', bbox_to_anchor=(1.01, 1.0), handletextpad=0.4)
axr.add_artist(leg)
axr.legend(handles=real_keys, fontsize=11, frameon=True, loc='lower right')
axr.set_xlabel(r'$X_0^{\rm mean}$', fontsize=20)
axr.set_ylabel(r'$F_\nu~\rm [mJy]$', fontsize=20)
_ticks(axr)
plt.savefig('plots/test_fluxflux_embedded_curves.png', dpi=200,
            bbox_inches='tight', bbox_extra_artists=(leg,))
plt.close()
print("Saved plots/test_fluxflux_embedded_curves.png")

# --- model versus the observed PG 0844 decomposition, mJy ------------------
if OBS is not None:
    print("\nModel (fiducial embedded mass) vs observed PG 0844, mJy. Model "
          "values are\ninterpolated in log-lambda onto the observed bands.")
    print(f"{'lam [A]':>9}{'host obs':>10}{'host mod':>10}{'ratio':>7}"
          f"{'bright obs':>12}{'bright mod':>12}{'faint obs':>11}"
          f"{'faint mod':>11}")

    def interp_band(vals):
        return np.exp(np.interp(np.log(OBS['lam']), np.log(bands),
                                np.log(np.maximum(vals, 1e-12))))

    h_mod = interp_band(rec_p)
    b_mod = interp_band(Fh_p[i_hi])
    f_mod = interp_band(Fh_p[i_lo])
    for k, lam in enumerate(OBS['lam']):
        print(f"{lam:>9.0f}{OBS['host'][k]:>10.3f}{h_mod[k]:>10.3f}"
              f"{h_mod[k]/max(OBS['host'][k], 1e-6):>7.2f}"
              f"{OBS['bright'][k]:>12.3f}{b_mod[k]:>12.3f}"
              f"{OBS['faint'][k]:>11.3f}{f_mod[k]:>11.3f}")

# --- Figure 4: line-of-sight view of the surface temperature ---------------
# Same rendering as transient/plot_los_view.py: the surface is projected onto
# the sky plane for the observer direction e = (sin i, 0, cos i), depth sorted
# so raised pillars cover what they hide, and colored by surface temperature.
LOS_LIM = 20.0                 # sky half-width, light days
LOS_VMIN, LOS_VMAX = 5.0e2, 4.0e4       # floor set by T_visc at the outer edge
_, _, tpk_fid = spot_t4(gp, M_EMBEDDED_FID)
T_los_spot = (T4_p[i_hi] + t4f) ** 0.25
T_los_nospot = T4_p[i_hi] ** 0.25
print(f"LOS map temperature range: no spots {T_los_nospot.min():.0f}-"
      f"{T_los_nospot.max():.0f} K, with spots {T_los_spot.min():.0f}-"
      f"{T_los_spot.max():.0f} K")
m_phys = M_EMBEDDED_FID * FLUX_SCALE
GEOM_NOTE = (rf'${len(disk.pillars)}~{{\rm pillars}},~h_p = {H_P}~{{\rm ld}},~'
             rf'\sigma_r = {SIGMA_R}~{{\rm ld}},~\sigma_\phi = {SIGMA_PHI}$')
_los_map(gp, T_los_spot, 'plots/test_fluxflux_embedded_tmap.png',
         notes=[GEOM_NOTE,
                r'$\rm embedded~stars:$ ' + mass_label(m_phys)
                + r' $\rm per~pillar,$ '
                + rf'$T_{{\rm spot}}^{{\rm peak}} = '
                  rf'{tpk_fid.min()/1e3:.1f}-{tpk_fid.max()/1e3:.1f}~{{\rm kK}}$'],
         lim=LOS_LIM, vmin=LOS_VMIN, vmax=LOS_VMAX)
_los_map(gp, T_los_nospot, 'plots/test_fluxflux_embedded_tmap_nospot.png',
         notes=[GEOM_NOTE, r'$\rm no~embedded~stars~(shadows~only)$'],
         lim=LOS_LIM, vmin=LOS_VMIN, vmax=LOS_VMAX)
