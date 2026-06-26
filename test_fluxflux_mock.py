"""Mock flux-flux decomposition from the pillar disc temperature profile.

Idea (Keith): the lamp varies, so each disc cell splits into a constant
viscous floor B_nu(T_visc(r)) plus a lamp-driven variable excess that
exists ONLY on illuminated cells. Shadowed cells (S=0) never see the
varying lamp, so their viscous emission lands entirely in the constant
("host-like") component of a flux-flux analysis. This script computes:

  F_total   : full SED with shadows (lamp on)            [black]
  F_const   : viscous-only SED, whole disc (lamp off)    [red]   = constant
  F_var     : F_total - F_const                          [blue]  = variable
  F_shadow  : viscous SED of the shadowed area only       [red --] = the
              red excess that masquerades as extra host light.

All integrals reuse the exact geometry/units of PillarDisk.compute_sed
(validated below), just with an explicit temperature map and an optional
per-cell weight so we can isolate the shadowed area.
"""
import numpy as np
import matplotlib.pyplot as plt

from pillardisk.pillar_disk import C, DAY, FNU_TO_MJY, H, K, ANGSTROM, PC_TO_LD
# Reuse the continuum-section disc + pillar builders (config_line.yaml,
# NGC 5548 fiducial). Importing is safe: the scan is under __main__.
from test_lag_spectrum import build_disk, add_pillars

# Fiducial continuum-section pillars: h_p=0.5, sigma_r=0.5, sigma_phi=0.3
# Thin pillars to suppress the optical bump: H_P=0.5, sigma_r=0.25, sigma_phi=0.3
# (the width-scan config giving a ~1.3x optical bump over the bowl).
H_P, SIGMA_R, SIGMA_PHI = 0.5, 0.25, 0.3
HLAMP = 50 * 0.00399          # lamp height = 50 r_g (r_g=0.00399 ld for 7e7 Msun)

# Place the disc at PG 0844+349's distance (z = 0.064, D_L ~ 287 Mpc for a
# flat H0=70, Om=0.3 cosmology) so the SED comes out in physical mJy. The
# model applies no (1+z) factor, so this is just the luminosity distance.
Z_PG0844, DL_MPC = 0.064, 287.0

wavelengths = np.logspace(np.log10(1000.0), np.log10(10000.0), 50)

disk = build_disk(hlamp=HLAMP)
add_pillars(disk, 100, H_P, SIGMA_R, SIGMA_PHI, r_min=5.0)    # pillar belt r>=5 ld
disk.d = DL_MPC * PC_TO_LD                 # override distance to PG 0844
print(f"Disc: rin={disk.rin:.3f} rout={disk.rout} ld, nr={disk.nr} "
      f"nphi={disk.nphi}, {len(disk.pillars)} pillars; "
      f"z={Z_PG0844}, D_L={DL_MPC} Mpc")

# --- Geometry and masks (on the disc's own grid) ---
r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')
h_2d = disk.get_height(r_2d, phi_2d)
S = disk._compute_shadow_mask(r_2d, phi_2d, h_2d)        # 1=lit, 0=shadow
f_shadow = 1.0 - S.mean(axis=1)                           # covering fraction(r)

# Viscous temperature map (axisymmetric T_visc(r), broadcast over phi)
tv_2d = np.interp(r_2d.ravel(), disk.r, disk.tv_base).reshape(r_2d.shape)
T_full = disk.get_temperature(r_2d, phi_2d, compute_shadows=True)


def _geom_weight():
    """Per-cell geometric weight da*max(dot,0) matching compute_sed,
    evaluated at the lower radius of each radial slab (shape nr-1, nphi)."""
    rl = disk.r[:-1]
    dr = np.diff(disk.r)
    dh = h_2d[1:, :] - h_2d[:-1, :]
    ds = np.sqrt(dr[:, None]**2 + dh**2)
    sintilt = dh / (ds + 1e-10)
    costilt = dr[:, None] / (ds + 1e-10)
    px = -sintilt * np.cos(disk.phi)[None, :]
    pz = costilt
    dot = disk.ex * px + disk.ez * pz
    dot = np.where(dot > 0, dot, 0.0)
    da = ds * rl[:, None] * disk.dphi
    return da * dot                                       # (nr-1, nphi)


WGEO = _geom_weight()
LD2CM = C * DAY
D_CM = disk.d * LD2CM
SED_UNITS = LD2CM**2 / D_CM**2 * FNU_TO_MJY    # solid-angle x mJy conversion


def b_nu(lam_A, T):
    """Physical Planck function B_nu in erg/cm^2/s/Hz/ster."""
    nu = C / (lam_A * ANGSTROM)
    x = np.clip(H * nu / (K * np.maximum(T, 1.0)), 1e-10, 700.0)
    return 2.0 * H * nu**3 / C**2 / (np.exp(x) - 1.0)


def sed_from_T(T_2d, weight=None, wl=None):
    """Physical F_nu (mJy) for an explicit temperature map, optional per-cell
    weight. Uses the validated compute_sed geometry (WGEO) but a true B_nu
    so the absolute scale is physical (compute_sed itself is F_lambda-like)."""
    wl = wavelengths if wl is None else np.atleast_1d(wl)
    Tlow = T_2d[:-1, :]
    w = WGEO if weight is None else WGEO * weight[:-1, :]
    f4 = disk.fcol**4
    out = np.empty(len(wl))
    for i, lam in enumerate(wl):
        B = b_nu(lam, Tlow * disk.fcol) / f4
        out[i] = np.sum(B * w)
    return out * SED_UNITS


# === commented out to speed up (remove the _SKIP guards to re-enable) ===
_SKIP_geo = r"""
# --- Geometry check: WGEO + units reproduce compute_sed's own (F_lambda-like)
# normalization to machine precision, confirming the integration is faithful. ---
chk_w = np.array([1500.0, 5100.0, 9000.0])
ref = disk.compute_sed(chk_w)
mine = np.array([np.sum(disk.planck_function(lam, T_full[:-1, :] * disk.fcol)
                        / disk.fcol**4 * WGEO) * SED_UNITS for lam in chk_w])
print("Geometry check vs compute_sed (its F_lambda-like units):")
for lam, a, b in zip(chk_w, ref, mine):
    print(f"  {lam:6.0f} A : compute_sed={a:.4e}  helper={b:.4e}  "
          f"rel diff={abs(a-b)/a:.2e}")
"""

# === commented out to speed up ===
_SKIP_dec = r"""
# --- Decomposition ---
F_total = sed_from_T(T_full)
F_const = sed_from_T(tv_2d)                  # viscous, whole disc (lamp off)
F_var = F_total - F_const                    # lamp-driven variable excess
F_shadow = sed_from_T(tv_2d, weight=1.0 - S)  # shadowed-area viscous light

# Report at V band (5100 A)
iV = np.argmin(np.abs(wavelengths - 5100.0))
print(f"\nAt 5100 A:  F_total={F_total[iV]:.3e}  F_var={F_var[iV]:.3e}  "
      f"F_const={F_const[iV]:.3e}  F_shadow={F_shadow[iV]:.3e} mJy")
print(f"  shadow / variable = {F_shadow[iV]/F_var[iV]:.2f}")
print(f"  shadow / total    = {F_shadow[iV]/F_total[iV]:.2f}")
print(f"  shadow covering fraction: {f_shadow.min():.3f}-{f_shadow.max():.3f}")
"""

# --- Plot ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def _ticks(ax):
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


# === commented out to speed up ===
_SKIP_fig1 = r"""
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Panel 1: SED decomposition
ax1.plot(wavelengths, F_total, color='black', lw=3, alpha=0.9,
         label=r'$\rm Total~(lamp~on)$')
ax1.plot(wavelengths, F_var, color='royalblue', lw=3, alpha=0.85,
         label=r'$\rm Variable~(irradiated)$')
ax1.plot(wavelengths, F_const, color='orangered', lw=3, alpha=0.85,
         label=r'$\rm Constant~(viscous,~all)$')
ax1.plot(wavelengths, F_shadow, color='orangered', lw=3, ls='--', alpha=0.85,
         label=r'$\rm Shadow~excess$')
# nu^{1/3} guide (F_nu propto lambda^{-1/3}) anchored to the variable curve
guide = F_var[iV] * (wavelengths / 5100.0)**(-1.0/3.0)
ax1.plot(wavelengths, guide, color='gray', lw=1.5, ls=':',
         label=r'$f_\nu \propto \nu^{1/3}$')
ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_xlabel(r'$\rm Wavelength~[\AA]$', fontsize=18)
ax1.set_ylabel(r'$F_\nu~\rm [mJy]$', fontsize=18)
ax1.set_ylim(top=F_total.max()*2, bottom=F_total.max()*1e-3)
ax1.legend(fontsize=13, frameon=False)
_ticks(ax1)

# Panel 2: inputs T_visc(r) and shadow covering fraction f_shadow(r)
ax2.plot(disk.r, disk.tv_base, color='orangered', lw=3, alpha=0.9)
ax2.set_xlabel(r'$r~\rm [light~days]$', fontsize=18)
ax2.set_ylabel(r'$T_{\rm visc}(r)~\rm [K]$', fontsize=18, color='orangered')
ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlim(disk.rin, disk.rout)
ax2.tick_params(axis='y', colors='orangered')
_ticks(ax2)
ax2b = ax2.twinx()
ax2b.plot(disk.r, f_shadow, color='royalblue', lw=3, alpha=0.9)
ax2b.set_ylabel(r'$f_{\rm shadow}(r)$', fontsize=18, color='royalblue')
ax2b.tick_params(axis='y', colors='royalblue', which='both', direction='in')
ax2b.tick_params(which='major', length=8, width=1.5, labelsize=14)
ax2b.tick_params(which='minor', length=4, width=1.0)
ax2b.minorticks_on()
ax2b.set_xlim(disk.rin, disk.rout)

plt.tight_layout()
plt.savefig('plots/test_fluxflux_mock.png', dpi=200, bbox_inches='tight')
print("\nSaved plots/test_fluxflux_mock.png")
"""

# ---------------------------------------------------------------------------
# Flux-flux diagram: vary the lamp luminosity and trace each band vs the UV
# driving band. T_irr scales as (lamp)^{1/4}, so scaling tx_base by `s`
# scales the lamp luminosity by s^4. Faint limit s->0 is the viscous (con-
# stant) floor; the band-vs-band loci extrapolate to that constant SED.
# ---------------------------------------------------------------------------
# PG 0844-like band set, UVW2 -> z (rest-frame effective wavelengths, A)
bands = np.array([1928.0, 2246.0, 2600.0, 3465.0, 4392.0,
                  5468.0, 6215.0, 7545.0, 8700.0])
# Diffuse floor: a constant (lamp-independent), wavelength-flat glow filling the
# pillar shadows, carrying F_DIFF of the local direct-irradiation power, with a
# Balmer bound-free rise blueward of 3646 A. SIGMA_I*T^4 is the bolometric
# blackbody intensity, spread flat in F_nu over the band frequency range DNU.
F_DIFF = 0.3
SIGMA_I = 2.0 * np.pi**4 * K**4 / (15.0 * H**3 * C**2)
DNU = C / (bands.min() * ANGSTROM) - C / (bands.max() * ANGSTROM)
band_lab = [r'$\rm UVW2$', r'$\rm UVM2$', r'$\rm UVW1$', r'$U$', r'$B$',
            r'$V$', r'$r$', r'$i$', r'$z$']
# wavelength colormap matching the PG 0844 plots (turbo, blue UV -> red IR)
band_col = list(plt.cm.turbo(np.linspace(0.06, 0.95, len(bands))))
iV = int(np.argmin(np.abs(bands - 5100.0)))          # V-ish band for scaling
# Lamp scale (T_irr ~ lamp^{1/4}). The lamp's full amplitude is [0, S_FULL],
# but a few-year campaign only catches it varying over a fraction CAMP_FRAC of
# that, the window [S_LO, S_HI]. Those campaign states define the X
# normalization and the linear fit; the fainter states (down to lamp off) are
# unobserved and only used for the true-host / anchor reference.
S_FULL = 1.35                                  # brightest lamp (full amplitude)
CAMP_FRAC = 0.20                               # fraction of the amplitude sampled
S_HI = S_FULL
S_LO = S_FULL * (1.0 - CAMP_FRAC)
s_vals = np.concatenate([
    np.linspace(0.0, S_LO, 8, endpoint=False),  # unobserved faint tail
    np.linspace(S_LO, S_HI, 14)])               # observed campaign window
CAMP_MIN = S_LO


def diffuse_sed(dk, wl):
    """Constant, wavelength-flat diffuse glow (mJy, no FLUX_SCALE) filling the
    pillar shadows: F_DIFF of the local direct-irradiation power, spread flat in
    F_nu with a Balmer rise blueward of 3646 A. Evaluated at the disc's current
    (nominal) tx_base, so it is the same in every lamp state. Zero for the bowl
    (no pillars -> no shadows)."""
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    tv2 = np.interp(r2.ravel(), dk.r, dk.tv_base).reshape(r2.shape)
    Tns = dk.get_temperature(r2, p2, compute_shadows=False)
    Tdir_nom = np.maximum(Tns**4 - tv2**4, 0.0)**0.25
    if len(dk.pillars) > 0:
        smask = dk._compute_shadow_mask(r2, p2, h2)
    else:
        smask = np.ones_like(r2)
    dh = h2[1:, :] - h2[:-1, :]
    ds = np.sqrt(np.diff(dk.r)[:, None]**2 + dh**2)
    px = -(dh / (ds + 1e-10)) * np.cos(dk.phi)[None, :]
    pz = np.diff(dk.r)[:, None] / (ds + 1e-10)
    dot = np.where(dk.ex * px + dk.ez * pz > 0, dk.ex * px + dk.ez * pz, 0.0)
    wgeo = ds * dk.r[:-1, None] * dk.dphi * dot
    units = (C*DAY)**2 / (dk.d*C*DAY)**2 * FNU_TO_MJY
    diff_level = F_DIFF * SIGMA_I * Tdir_nom[:-1, :]**4 / DNU
    D0 = np.sum(diff_level * (1.0 - smask[:-1, :]) * wgeo)      # flat-F_nu amplitude
    wl = np.atleast_1d(wl)
    balmer = np.where(wl < 3646.0, 1.5, 1.0)
    return D0 * balmer * units


def fluxflux(dk, bands, s_vals):
    """Sweep the lamp and return per-band band fluxes (n_s, n_band) plus the
    faint-limit constant SED, for an arbitrary disc (pillared or bowl)."""
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    rl = dk.r[:-1]
    dh = h2[1:, :] - h2[:-1, :]
    ds = np.sqrt(np.diff(dk.r)[:, None]**2 + dh**2)
    px = -(dh / (ds + 1e-10)) * np.cos(dk.phi)[None, :]
    pz = np.diff(dk.r)[:, None] / (ds + 1e-10)
    dot = np.where(dk.ex * px + dk.ez * pz > 0, dk.ex * px + dk.ez * pz, 0.0)
    wgeo = ds * rl[:, None] * dk.dphi * dot
    units = (C*DAY)**2 / (dk.d*C*DAY)**2 * FNU_TO_MJY

    def sed(T2d, wl):
        Tl = T2d[:-1, :]
        return np.array([np.sum(b_nu(l, Tl*dk.fcol) / dk.fcol**4 * wgeo)
                         for l in wl]) * units

    diffuse_band = diffuse_sed(dk, bands)          # constant diffuse glow (nominal state)
    tv2 = np.interp(r2.ravel(), dk.r, dk.tv_base).reshape(r2.shape)
    tx0 = dk.tx_base.copy()
    F = np.zeros((len(s_vals), len(bands)))
    for j, s in enumerate(s_vals):
        dk.tx_base = tx0 * s
        F[j] = sed(dk.get_temperature(r2, p2, compute_shadows=True), bands) + diffuse_band
    dk.tx_base = tx0
    return F, sed(tv2, bands) + diffuse_band


def eta_lp(dk):
    """Bolometric lamp/viscous power ratio at the nominal (s=1) state:
    eta_LP = Int(T_irr^4 dA) / Int(T_visc^4 dA) over the disc surface."""
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    ds = np.sqrt(np.diff(dk.r)[:, None]**2 + (h2[1:, :] - h2[:-1, :])**2)
    da = ds * dk.r[:-1, None] * dk.dphi                 # surface area element
    tv2 = np.interp(r2.ravel(), dk.r, dk.tv_base).reshape(r2.shape)
    T = dk.get_temperature(r2, p2, compute_shadows=True)
    tx4 = np.maximum(T**4 - tv2**4, 0.0)                 # irradiation T_irr^4
    return np.sum(tx4[:-1, :] * da) / np.sum(tv2[:-1, :]**4 * da)


# Bowl baseline: same disc, no pillars, same PG 0844 distance
bowl = build_disk(hlamp=HLAMP)
bowl.d = DL_MPC * PC_TO_LD

# Rescale the irradiation so the bolometric lamp power = ETA_LP x the viscous
# power at the nominal (s=1) state. config_line.yaml has tirrad_tvisc_ratio=10,
# i.e. eta_LP ~ 10^4 (viscous negligible); eta_LP=3 is far more physical and
# makes the viscous SED a meaningful fraction of the optical flux.
ETA_LP = 50.0
eta0 = eta_lp(disk)
f_eta = (ETA_LP / eta0) ** 0.25
disk.tx_base *= f_eta
bowl.tx_base *= f_eta
print(f"eta_LP: original {eta0:.1f} -> rescaled to {eta_lp(disk):.2f} "
      f"(factor {f_eta:.3f} on T_irr)")
disk.tv_base *= 1.5            # grow the viscous floor (lamp/viscous ratio 50 -> ~10)
bowl.tv_base *= 1.5
print(f"viscous floor raised x1.5 -> lamp/viscous ratio now {eta_lp(disk):.1f}")

F_states, F_const_b = fluxflux(disk, bands, s_vals)   # pillared disc
F_bowl, F_const_bowl = fluxflux(bowl, bands, s_vals)

# The absolute SED scale is arbitrary; rescale so the bright-state optical AGN
# lands at a realistic ~few-mJy level (PG 0844 bright-state V AGN ~ 9 mJy).
FLUX_SCALE = 9.0 / F_states[:, iV].max()
F_states *= FLUX_SCALE
F_bowl *= FLUX_SCALE
F_const_b *= FLUX_SCALE
F_const_bowl *= FLUX_SCALE

# Flux-flux x-axis: the driving light curve built the PyROA way. PyROA fits
# each band as F_i = A_i X(t) + B_i and merges them, so the driver is the mean
# of ALL bands after standardizing each (subtract mean B_i, divide by rms A_i)
# over the CAMPAIGN states -- NOT the shortest-wavelength band. The full sweep
# (incl. the non-variable tail) is mapped through the campaign statistics, so
# the loci extend to the floor while the linear fit is anchored to the campaign.
camp = s_vals >= CAMP_MIN


def driving_lc(F):
    """PyROA-style driving light curve: standardize each band over the
    campaign, average across all bands, then renormalize to zero mean / unit
    variance over the campaign."""
    B = F[camp].mean(axis=0)
    A = F[camp].std(axis=0)
    X = ((F - B) / A).mean(axis=1)
    return (X - X[camp].mean()) / X[camp].std()


Xdrive = driving_lc(F_states)
Xbowl = driving_lc(F_bowl)
X_const = Xdrive[0]            # lamp off (s=0) -> non-variable point
X_const_bowl = Xbowl[0]


def nonvar_anchor(fit):
    """X_non-var: the least-negative band F=0 crossing -- the FIRST band to
    reach zero flux as the driver decreases (the host-negligible blue band,
    as anchored on UVW2 in the PG 0844 flux-flux). Evaluating every band here
    gives a non-negative recovered constant by construction."""
    return np.max(-fit[:, 1] / fit[:, 0])


def draw_ff(ax):
    """Flux-flux loci (pillars solid, bowl dashed) extended to the non-variable
    end, campaign linear fits, the recovered constant SED at X_non-var
    (star=pillars, circle=bowl), and the true viscous floor (squares)."""
    pfit = np.array([np.polyfit(Xdrive[camp], F_states[camp, k], 1) for k in range(len(bands))])
    bfit = np.array([np.polyfit(Xbowl[camp], F_bowl[camp, k], 1) for k in range(len(bands))])
    Xnv_p = nonvar_anchor(pfit)            # pillar non-variability anchor
    Xnv_b = nonvar_anchor(bfit)            # bowl non-variability anchor
    for k in range(len(bands)):
        c = band_col[k]
        # full loci, extended down to the non-variable end
        ax.plot(Xdrive, F_states[:, k], 'o-', color=c, lw=2.5, ms=6,
                alpha=0.9, label=band_lab[k])
        ax.plot(Xbowl, F_bowl[:, k], '--', color=c, lw=2.0, alpha=0.7)
        # campaign linear fits, extrapolated to X_non-var
        xp = np.array([Xnv_p, Xdrive.max()])
        ax.plot(xp, pfit[k, 0] * xp + pfit[k, 1], ':', color=c, lw=1.2, alpha=0.6)
        xb = np.array([Xnv_b, Xbowl.max()])
        ax.plot(xb, bfit[k, 0] * xb + bfit[k, 1], '-.', color=c, lw=1.0, alpha=0.5)
        # recovered constant at X_non-var: star=pillars, circle=bowl
        ax.plot(Xnv_p, pfit[k, 0] * Xnv_p + pfit[k, 1], marker='*', color=c,
                ms=17, mec='k', mew=0.9, ls='none', zorder=5)
        ax.plot(Xnv_b, bfit[k, 0] * Xnv_b + bfit[k, 1], marker='o',
                mfc='white', mec=c, mew=2.0, ms=10, ls='none', zorder=5)
        # true viscous floor = left end of the extended loci (lamp off)
        ax.plot(X_const, F_const_b[k], marker='s', color=c, ms=8,
                mec='k', mew=0.9, ls='none', zorder=6)
    ax.axhline(0.0, color='gray', ls='--', lw=1.5, alpha=0.7)
    ax.axvline(Xnv_p, color='gray', ls='--', lw=1.5, alpha=0.7)
    ax.axvline(Xnv_b, color='gray', ls=':', lw=1.5, alpha=0.7)
    return Xnv_p, Xnv_b


from matplotlib.lines import Line2D

# === commented out to speed up ===
_SKIP_curves_driver = r"""
fig2, (ax, axz) = plt.subplots(1, 2, figsize=(15, 7))
X0p, X0b = draw_ff(ax)
draw_ff(axz)
X0 = min(X0p, X0b)

# Full panel: legends + labels
style_keys = [Line2D([], [], color='k', ls='-', marker='o', label=r'$\rm pillars$'),
              Line2D([], [], color='k', ls='--', label=r'$\rm bowl$'),
              Line2D([], [], color='k', marker='*', ms=14, mec='k', ls='none',
                     label=r'$\rm const~(pillar)$'),
              Line2D([], [], color='k', marker='o', mfc='white', mec='k', ms=9,
                     ls='none', label=r'$\rm const~(bowl)$'),
              Line2D([], [], color='k', marker='s', mec='k', ms=8, ls='none',
                     label=r'$\rm true~floor$')]
leg1 = ax.legend(fontsize=14, frameon=False, title=r'$\rm band$', loc='upper left')
ax.add_artist(leg1)
ax.legend(handles=style_keys, fontsize=13, frameon=False, loc='lower right')
ax.set_xlabel(r'$X_0^{\rm mean}$', fontsize=20)
ax.set_ylabel(r'$F_\nu~\rm [mJy]$', fontsize=20)
_ticks(ax)

# Zoom panel: lower-left, around the intercepts and the constant stars
axz.set_xlim(X0 - 0.10, -0.3)
axz.set_ylim(-2.0, 0.25 * F_states.max())
axz.set_xlabel(r'$X_0^{\rm mean}$', fontsize=20)
axz.set_ylabel(r'$F_\nu~\rm [mJy]$', fontsize=20)
_ticks(axz)

plt.tight_layout()
plt.savefig('plots/test_fluxflux_curves.png', dpi=200, bbox_inches='tight')
print("Saved plots/test_fluxflux_curves.png")
pfit = np.array([np.polyfit(Xdrive[camp], F_states[camp, k], 1) for k in range(len(bands))])
bfit = np.array([np.polyfit(Xbowl[camp], F_bowl[camp, k], 1) for k in range(len(bands))])
Xnv_p, Xnv_b = nonvar_anchor(pfit), nonvar_anchor(bfit)
print(f"\nPyROA-style driving LC (mean of all bands): X range "
      f"[{Xdrive.min():.2f}, {Xdrive.max():.2f}]; X_non-var(pillar) = {Xnv_p:.2f}")
print("Recovered constant SED at X_non-var (first band to reach F=0), "
      "vs the true viscous floor (mJy):")
for k, lam in enumerate(bands):
    rec_p = pfit[k, 0]*Xnv_p + pfit[k, 1]
    rec_b = bfit[k, 0]*Xnv_b + bfit[k, 1]
    print(f"  {lam:6.0f} A : recovered(pillar)={rec_p:7.3f}  "
          f"recovered(bowl)={rec_b:7.3f}  true_floor={F_const_b[k]:.3f}")

# ---------------------------------------------------------------------------
# Driver-choice comparison: same loci, but with X built from the underlying
# DRIVING signal (the lamp) rather than the 1500 A band, as in Keith's
# bbdisk reference. Two separate figures: L_lamp (~ s^4) and T_irr (~ s).
# Pillar and bowl share the same driver X (the lamp is disc-independent).
# X is normalized over the full sweep so the lamp-off floor lands near X~-1.
# ---------------------------------------------------------------------------
style_keys2 = [Line2D([], [], color='k', ls='-', marker='o', label=r'$\rm pillars$'),
               Line2D([], [], color='k', ls='--', label=r'$\rm bowl$'),
               Line2D([], [], color='k', marker='*', mec='k', ms=14, ls='none',
                      label=r'$\rm recovered$'),
               Line2D([], [], color='k', marker='s', mec='k', ms=8, ls='none',
                      label=r'$\rm true~floor$')]


def make_driver_fig(D, fname, xlabel):
    X = (D - D.mean()) / D.std()
    Xfloor = X[0]                     # s=0 lamp off = non-variable point
    fig, axx = plt.subplots(figsize=(8, 7))
    for k in range(len(bands)):
        c = band_col[k]
        axx.plot(X, F_states[:, k], 'o-', color=c, lw=2.5, ms=6, alpha=0.9,
                 label=band_lab[k])
        axx.plot(X, F_bowl[:, k], '--', color=c, lw=2.0, alpha=0.65)
        m, b = np.polyfit(X[camp], F_states[camp, k], 1)
        xext = np.array([Xfloor, X.max()])
        axx.plot(xext, m * xext + b, ':', color=c, lw=1.1, alpha=0.55)
        axx.plot(Xfloor, m * Xfloor + b, marker='*', color=c, ms=16, mec='k',
                 mew=0.8, ls='none', zorder=5)
        axx.plot(Xfloor, F_const_b[k], marker='s', color=c, ms=8, mec='k',
                 mew=0.8, ls='none', zorder=6)
    axx.axhline(0.0, color='gray', ls='--', lw=1.5, alpha=0.7)
    axx.axvline(Xfloor, color='gray', ls='--', lw=1.5, alpha=0.7)
    leg = axx.legend(fontsize=13, frameon=False, title=r'$\rm band$', loc='upper left')
    axx.add_artist(leg)
    axx.legend(handles=style_keys2, fontsize=12, frameon=False, loc='lower right')
    axx.set_xlabel(xlabel, fontsize=20)
    axx.set_ylabel(r'$F_\nu~\rm [mJy]$', fontsize=20)
    _ticks(axx)
    plt.tight_layout()
    plt.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {fname}")


make_driver_fig(s_vals**4, 'plots/test_fluxflux_driver_L.png', r'$X_0^{L}$')
make_driver_fig(s_vals, 'plots/test_fluxflux_driver_T.png', r'$X_0^{T}$')
"""

# ---------------------------------------------------------------------------
# Add a constant host SED to every band, re-derive the PyROA driving light
# curve, and plot separately. Because the driver standardizes each band
# (subtract B_i), a constant host cancels and X_0 is unchanged -- the point is
# to show that the flux-flux extrapolation then recovers ~the host.
# Host model: a generic red power-law galaxy SED F_nu ~ nu^{-2} (~lambda^2),
# normalized to HOST_5100 mJy at 5100 A (old-population host, not the real
# PG 0844 host).
# beta=3 (red host), HOST_5100=1.825 -> host fraction ~0.2 at V (5468 A).
HOST_BETA, HOST_5100 = 3.0, 1.825
HOST = HOST_5100 * (bands / 5100.0) ** HOST_BETA
F_states_h = F_states + HOST
F_bowl_h = F_bowl + HOST
Xh = driving_lc(F_states_h)            # identical to Xdrive: host cancels
Xhb = driving_lc(F_bowl_h)
F_const_h = F_const_b + HOST           # true constant = host + viscous floor

# === commented out to speed up ===
_SKIP_host = r"""
fig4, axh = plt.subplots(figsize=(8, 7))
ph = np.array([np.polyfit(Xh[camp], F_states_h[camp, k], 1) for k in range(len(bands))])
bh = np.array([np.polyfit(Xhb[camp], F_bowl_h[camp, k], 1) for k in range(len(bands))])
# X_non-var = first band to reach F=0 (anchors the host to a non-negative SED)
Xh_nv, Xhb_nv = nonvar_anchor(ph), nonvar_anchor(bh)
for k in range(len(bands)):
    c = band_col[k]
    axh.plot(Xh, F_states_h[:, k], 'o-', color=c, lw=2.5, ms=6, alpha=0.9,
             label=band_lab[k])
    axh.plot(Xhb, F_bowl_h[:, k], '--', color=c, lw=2.0, alpha=0.65)
    xp = np.array([Xh_nv, Xh.max()])
    axh.plot(xp, ph[k, 0] * xp + ph[k, 1], ':', color=c, lw=1.2, alpha=0.6)
    axh.plot(Xh_nv, ph[k, 0] * Xh_nv + ph[k, 1], marker='*', color=c,
             ms=16, mec='k', mew=0.8, ls='none', zorder=5)
    axh.plot(Xh_nv, F_const_h[k], marker='s', color=c, ms=8, mec='k',
             mew=0.8, ls='none', zorder=6)
axh.axhline(0.0, color='gray', ls='--', lw=1.5, alpha=0.7)
axh.axvline(Xh_nv, color='gray', ls='--', lw=1.5, alpha=0.7)
host_keys = [Line2D([], [], color='k', ls='-', marker='o', label=r'$\rm pillars$'),
             Line2D([], [], color='k', ls='--', label=r'$\rm bowl$'),
             Line2D([], [], color='k', marker='*', mec='k', ms=14, ls='none',
                    label=r'$\rm recovered$'),
             Line2D([], [], color='k', marker='s', mec='k', ms=8, ls='none',
                    label=r'$\rm true~host$')]
leg = axh.legend(fontsize=13, frameon=False, title=r'$\rm band$', loc='upper left')
axh.add_artist(leg)
axh.legend(handles=host_keys, fontsize=12, frameon=False, loc='lower right')
axh.set_xlabel(r'$X_0^{\rm mean}$', fontsize=20)
axh.set_ylabel(r'$F_\nu~\rm [mJy]$', fontsize=20)
_ticks(axh)
plt.tight_layout()
plt.savefig('plots/test_fluxflux_host.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved plots/test_fluxflux_host.png")
print(f"Host added -> recovered host at X_non-var={Xh_nv:.2f} vs true host (mJy):")
for k, lam in enumerate(bands):
    rec = ph[k, 0] * Xh_nv + ph[k, 1]
    print(f"  {lam:6.0f} A : recovered={rec:6.2f}  true_host={F_const_h[k]:.2f}  "
          f"(input host={HOST[k]:.2f})")
"""

# ---------------------------------------------------------------------------
# Realistic campaign: the faint / lamp-off states are never observed, so fit
# the linear flux-flux relation only over the observed campaign (camp, the same
# window the driver is normalized over), then extrapolate to recover the host.
# The unobserved fainter tail is drawn transparent. Recovered host shown for
# pillars (stars) and bowl (circles); true host as squares.
# ---------------------------------------------------------------------------
# Observed campaign = the normalization window (camp): the driving LC has zero
# mean / unit variance over it, so X_faint..X_bright are its actual edges
# (centered on 0), not an arbitrary cut. Fainter states (s < CAMP_MIN, down to
# lamp off) are unobserved and drawn transparent.
phr = np.array([np.polyfit(Xh[camp], F_states_h[camp, k], 1) for k in range(len(bands))])
bhr = np.array([np.polyfit(Xhb[camp], F_bowl_h[camp, k], 1) for k in range(len(bands))])
Xnv_pr, Xnv_br = nonvar_anchor(phr), nonvar_anchor(bhr)
Xf_p, Xb_p = Xh[camp].min(), Xh[camp].max()       # pillar observed edges
Xf_b, Xb_b = Xhb[camp].min(), Xhb[camp].max()     # bowl observed edges

fig5, axr = plt.subplots(figsize=(8, 7))
for k in range(len(bands)):
    c = band_col[k]
    # loci: unobserved faint tail transparent, observed campaign opaque
    axr.plot(Xh, F_states_h[:, k], 'o-', color=c, lw=2.5, ms=6, alpha=0.22)
    axr.plot(Xh[camp], F_states_h[camp, k], 'o-', color=c, lw=2.5, ms=6,
             alpha=0.9, label=band_lab[k])
    axr.plot(Xhb, F_bowl_h[:, k], marker='s', ls='--', color=c, mfc='white',
             lw=2.0, ms=5, alpha=0.18)
    axr.plot(Xhb[camp], F_bowl_h[camp, k], marker='s', ls='--', color=c,
             mfc='white', lw=2.0, ms=5, alpha=0.6)
    # pillar fit: opaque over the observed range, transparent in the extrapolation
    xin = np.array([Xf_p, Xb_p])
    axr.plot(xin, phr[k, 0]*xin + phr[k, 1], ':', color=c, lw=1.5, alpha=0.7)
    axr.plot([Xnv_pr, Xf_p], phr[k, 0]*np.array([Xnv_pr, Xf_p]) + phr[k, 1],
             ':', color=c, lw=1.2, alpha=0.3)
    # bowl fit (dash-dot), same opaque/transparent split
    xib = np.array([Xf_b, Xb_b])
    axr.plot(xib, bhr[k, 0]*xib + bhr[k, 1], '-.', color=c, lw=1.3, alpha=0.6)
    axr.plot([Xnv_br, Xf_b], bhr[k, 0]*np.array([Xnv_br, Xf_b]) + bhr[k, 1],
             '-.', color=c, lw=1.0, alpha=0.25)
    # recovered host: star=pillars, circle=bowl
    axr.plot(Xnv_pr, phr[k, 0]*Xnv_pr + phr[k, 1], marker='*', color=c,
             ms=16, mec='k', mew=0.8, ls='none', zorder=5)
    axr.plot(Xnv_br, bhr[k, 0]*Xnv_br + bhr[k, 1], marker='o', mfc='white',
             mec=c, mew=2.0, ms=10, ls='none', zorder=5)
    # true host (input SED, same for both) marked at the pillar anchor
    axr.plot(Xnv_pr, F_const_h[k], marker='D', color=c, ms=8, mec='k',
             mew=0.8, ls='none', zorder=6)
axr.axhline(0.0, color='gray', ls='--', lw=1.5, alpha=0.7)
axr.axvline(Xnv_pr, color='gray', ls='--', lw=1.5, alpha=0.7)     # pillar anchor
axr.axvline(Xnv_br, color='gray', ls='-.', lw=1.2, alpha=0.6)     # bowl anchor
# observed-range edges X_faint and X_bright (campaign min/max of the driver)
X_faint, X_bright = Xf_p, Xb_p
ytop = axr.get_ylim()[1]
for xv, lab in [(X_faint, r'$X_{\rm faint}$'), (X_bright, r'$X_{\rm bright}$')]:
    axr.axvline(xv, color='k', ls=':', lw=1.5, alpha=0.6)
    axr.text(xv, ytop, lab, rotation=90, va='top', ha='right',
             fontsize=11, color='k')
real_keys = [Line2D([], [], color='k', ls='-', marker='o', label=r'$\rm pillars$'),
             Line2D([], [], color='k', ls='--', marker='s', mfc='white',
                    label=r'$\rm bowl$'),
             Line2D([], [], color='k', marker='*', mec='k', ms=14, ls='none',
                    label=r'$\rm recovered,~pillar$'),
             Line2D([], [], color='k', marker='o', mfc='white', mec='k', ms=9,
                    ls='none', label=r'$\rm recovered,~bowl$'),
             Line2D([], [], color='k', marker='D', mec='k', ms=8, ls='none',
                    label=r'$\rm true~host$')]
# Band legend on the right, single column, reversed wavelength order (z -> UVW2)
band_handles = [Line2D([], [], color=band_col[k], marker='o', ls='-',
                       label=band_lab[k]) for k in range(len(bands))][::-1]
leg = axr.legend(handles=band_handles, fontsize=12, frameon=True, ncol=1,
                 loc='upper left', bbox_to_anchor=(1.01, 1.0),
                 handletextpad=0.4)
axr.add_artist(leg)
axr.legend(handles=real_keys, fontsize=11, frameon=True, loc='lower right')
axr.set_xlabel(r'$X_0^{\rm mean}$', fontsize=20)
axr.set_ylabel(r'$F_\nu~\rm [mJy]$', fontsize=20)
_ticks(axr)
plt.savefig('plots/test_fluxflux_realistic.png', dpi=200,
            bbox_inches='tight', bbox_extra_artists=(leg,))
plt.close()
print(f"Saved plots/test_fluxflux_realistic.png  (X_faint={X_faint:.2f}, "
      f"X_bright={X_bright:.2f}, X_non-var pillar={Xnv_pr:.2f} bowl={Xnv_br:.2f})")
print("Recovered host (pillar / bowl) vs true host (mJy):")
for k, lam in enumerate(bands):
    rp = phr[k, 0]*Xnv_pr + phr[k, 1]
    rb = bhr[k, 0]*Xnv_br + bhr[k, 1]
    print(f"  {lam:6.0f} A : pillar={rp:6.2f}  bowl={rb:6.2f}  true_host={F_const_h[k]:.2f}")

# ---------------------------------------------------------------------------
# SED view (F_nu vs wavelength) of the decomposition: the host galaxy SED, the
# fluxflux-recovered host from the pillar and bowl fits, the high/low (bright/
# faint) host-subtracted AGN states, and the true viscous SED.
# ---------------------------------------------------------------------------
i_hi = len(s_vals) - 1                    # brightest state
i_lo = int(np.argmax(camp))               # faintest observed (campaign) state
rec_p = phr[:, 0]*Xnv_pr + phr[:, 1]      # recovered host, pillar fit
rec_b = bhr[:, 0]*Xnv_br + bhr[:, 1]      # recovered host, bowl fit
sB, sF = s_vals[i_hi], s_vals[i_lo]
WL = np.logspace(np.log10(800), np.log10(20000), 80)    # AGN curves: 800 A -> 2 micron
nuWL = C / (WL * ANGSTROM)
WL_host = np.logspace(np.log10(1000), np.log10(20000), 70)   # host kept on its valid range
nuWL_host = C / (WL_host * ANGSTROM)
nub = C / (bands * ANGSTROM)


def sed_disc(dk, wl, s):
    """True AGN SED (mJy, host-free) for disc dk at lamp state s, scaled to the
    plot's flux normalization."""
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    dsg = np.sqrt(np.diff(dk.r)[:, None]**2 + (h2[1:, :] - h2[:-1, :])**2)
    px = -(h2[1:, :] - h2[:-1, :]) / (dsg + 1e-10) * np.cos(dk.phi)[None, :]
    pz = np.diff(dk.r)[:, None] / (dsg + 1e-10)
    dot = np.where(dk.ex*px + dk.ez*pz > 0, dk.ex*px + dk.ez*pz, 0.0)
    wg = dsg * dk.r[:-1, None] * dk.dphi * dot
    u = (C*DAY)**2 / (dk.d*C*DAY)**2 * FNU_TO_MJY * FLUX_SCALE
    t0 = dk.tx_base.copy(); dk.tx_base = t0 * s
    T = dk.get_temperature(r2, p2, compute_shadows=True); dk.tx_base = t0
    return np.array([np.sum(b_nu(l, T[:-1, :]*dk.fcol) / dk.fcol**4 * wg)
                     for l in wl]) * u + diffuse_sed(dk, wl) * FLUX_SCALE


TpB, TpF = sed_disc(disk, WL, sB), sed_disc(disk, WL, sF)   # pillar true
TbB, TbF = sed_disc(bowl, WL, sB), sed_disc(bowl, WL, sF)   # bowl true
HOST_WL = HOST_5100 * (WL_host / 5100.0) ** HOST_BETA      # host on its 1000-20000 A range
# inferred AGN at the bands = observed (AGN+host) - recovered host
inf_pB = F_states_h[i_hi] - rec_p
inf_pF = F_states_h[i_lo] - rec_p
inf_bB = F_bowl_h[i_hi] - rec_b
inf_bF = F_bowl_h[i_lo] - rec_b

# nu L_nu [erg s^-1] = 4 pi D^2 * nu F_nu. F_nu is in mJy (1e-26 erg s^-1 cm^-2 Hz^-1).
NUL_SCALE = 1e-26 * 4.0 * np.pi * D_CM**2

fig6, axs = plt.subplots(figsize=(8.5, 7))
# true AGN (larger range) and true host -- lines
axs.plot(WL, nuWL*TpB*NUL_SCALE, '-', color='royalblue', lw=2.5, label=r'$\rm high~AGN~(pillars)$')
axs.plot(WL, nuWL*TpF*NUL_SCALE, '-', color='navy', lw=2.5, label=r'$\rm low~AGN~(pillars)$')
axs.plot(WL, nuWL*TbB*NUL_SCALE, '--', color='royalblue', lw=2.0, alpha=0.6, label=r'$\rm high~AGN~(bowl)$')
axs.plot(WL, nuWL*TbF*NUL_SCALE, '--', color='navy', lw=2.0, alpha=0.6, label=r'$\rm low~AGN~(bowl)$')
axs.plot(WL_host, nuWL_host*HOST_WL*NUL_SCALE, '-', color='firebrick', lw=2.5, label=r'$\rm host~(true)$')
# recovered host at the bands (flux-flux fit)
mp = rec_p > 0
mb = rec_b > 0
axs.plot(bands[mp], nub[mp]*rec_p[mp]*NUL_SCALE, '*', color='orange', ms=16, mec='k',
         mew=0.8, ls='none', zorder=12, label=r'$\rm host~(pillar~fit)$')
axs.plot(bands[mb], nub[mb]*rec_b[mb]*NUL_SCALE, 'D', mfc='white', mec='darkorange',
         mew=2.0, ms=9, ls='none', label=r'$\rm host~(bowl~fit)$')
# inferred AGN at the bands -- markers (x = pillars, + = bowl; blue=high, navy=low)
axs.plot(bands, nub*inf_pB*NUL_SCALE, 'x', color='royalblue', ms=9, mew=2.2, ls=':', lw=1.3)
axs.plot(bands, nub*inf_pF*NUL_SCALE, 'x', color='navy', ms=9, mew=2.2, ls=':', lw=1.3)
axs.plot(bands, nub*inf_bB*NUL_SCALE, '+', color='royalblue', ms=11, mew=2.2, ls=':', lw=1.3)
axs.plot(bands, nub*inf_bF*NUL_SCALE, '+', color='navy', ms=11, mew=2.2, ls=':', lw=1.3)
inf_keys = [Line2D([], [], color='k', marker='x', ms=9, mew=2.2, ls='none',
                   label=r'$\rm inferred~AGN~(pillars)$'),
            Line2D([], [], color='k', marker='+', ms=11, mew=2.2, ls='none',
                   label=r'$\rm inferred~AGN~(bowl)$')]
axs.set_xscale('log'); axs.set_yscale('log')
allnu = np.concatenate([nuWL*TpB, nuWL*TbF, nuWL_host*HOST_WL]) * NUL_SCALE
axs.set_ylim(0.4 * allnu.min(), 3.0 * allnu.max())
axs.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=20)
axs.set_ylabel(r'$\nu L_\nu~\rm [erg~s^{-1}]$', fontsize=20)
leg1 = axs.legend(fontsize=11, frameon=True, loc='upper right', ncol=2)
axs.add_artist(leg1)
axs.legend(handles=inf_keys, fontsize=11, frameon=False, loc='upper left')
_ticks(axs)
plt.tight_layout()
plt.savefig('plots/test_fluxflux_sed.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved plots/test_fluxflux_sed.png")

# ---------------------------------------------------------------------------
# T(r) profile with pillars (bright state): azimuthal mean plus the 5-95th
# percentile spread over phi, which shows the cool shadow lanes (low) and the
# hot illuminated pillar walls (high) against the smooth bowl profile.
# ---------------------------------------------------------------------------
# === T(r) profile plot (enabled) ===
sB = s_vals[-1]
r2t, p2t = np.meshgrid(disk.r, disk.phi, indexing='ij')
tx0 = disk.tx_base.copy(); disk.tx_base = tx0 * sB
Tp = disk.get_temperature(r2t, p2t, compute_shadows=True); disk.tx_base = tx0
r2b, p2b = np.meshgrid(bowl.r, bowl.phi, indexing='ij')
tx0b = bowl.tx_base.copy(); bowl.tx_base = tx0b * sB
Tb = bowl.get_temperature(r2b, p2b, compute_shadows=True).mean(axis=1)
bowl.tx_base = tx0b
Tmean = Tp.mean(axis=1)
Tlo = np.percentile(Tp, 5, axis=1)
Thi = np.percentile(Tp, 95, axis=1)

fig7, axt = plt.subplots(figsize=(8, 7))
axt.fill_between(disk.r, Tlo, Thi, color='orange', alpha=0.3, lw=0,
                 label=r'$\rm pillars:~shadow{-}wall~range$')
axt.plot(disk.r, Tmean, color='orangered', lw=2.5, label=r'$\rm pillars~(mean)$')
axt.plot(bowl.r, Tb, color='navy', lw=2.5, ls='--', label=r'$\rm bowl$')
axt.set_xscale('linear'); axt.set_yscale('log')
axt.set_xlim(disk.rin, disk.rout)
axt.set_xlabel(r'$r~\rm [light~days]$', fontsize=20)
axt.set_ylabel(r'$T~\rm [K]$', fontsize=20)
axt.legend(fontsize=13, frameon=False, loc='lower left')
_ticks(axt)
plt.tight_layout()
plt.savefig('plots/test_fluxflux_tr.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved plots/test_fluxflux_tr.png")
