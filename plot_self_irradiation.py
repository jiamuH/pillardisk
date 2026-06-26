r"""Self-irradiation of the bowl by the hot inner disc, for several rim heights.

Two panels:
  Top: log Phi vs r -- lamp irradiation plus the inner-disc irradiation
       received by the bowl surface, one curve per rim height h1
       (bowl z = h1 (r/r0)^beta, beta = 10).
  Bottom: ratio F_disc / F_lamp vs r per h1.

Inner disc is treated as a Lambertian point source at the origin with
one-sided ionizing photon rate Q_disc = int Phi_exact(r) 2 pi r dr
(exact Planck photon integral above 1 Ryd, no-torque T profile):

  F_disc = Q_disc cos(theta_e) cos(theta_inc) / (pi d^2),  cos(theta_e) ~ z/r
  F_lamp = Q_lamp cos(theta_inc) / (4 pi d^2)
  F_disc / F_lamp = 4 (Q_disc / Q_lamp) cos(theta_e)

(The receiving-side cos(theta_inc) cancels in the ratio; the top panel
plots both fluxes for a face-on flat receiving element, cos(theta_inc)
identical for the two sources.)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
    'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

# Constants (cgs)
H_PL = 6.626e-27
K_B = 1.381e-16
C_CGS = 2.998e10
NU_0 = 13.6 * 1.602e-12 / H_PL
G_CGS = 6.674e-8
M_SUN = 1.989e33
LD_CM = C_CGS * 86400.0

# Parameters (config_line.yaml)
M_BH = 7.0e7
R_ISCO_RG = 2.0
HLAMP = 0.05
LOG_PHI_INNER = 25.0
TV1 = 446.0
R0 = 20.0
BETA = 100.0
R_OUT = 20.0
R_PILLAR_MIN = 2.0

H1_LIST = [0.05, 0.15, 0.5, 5.0]
COLORS = ['orangered', 'forestgreen', 'darkviolet', 'saddlebrown']

r_g_ld = G_CGS * M_BH * M_SUN / C_CGS**2 / LD_CM
r_in = R_ISCO_RG * r_g_ld

# ---- Lamp photon rate ----
d_in = np.sqrt(r_in**2 + HLAMP**2)
cos_in = HLAMP / d_in
Q_lamp = 10.0**LOG_PHI_INNER * 4.0 * np.pi * (d_in * LD_CM)**2 / cos_in

# ---- Inner-disc photon rate (exact Planck integral, no-torque T) ----
def wien_integral(x0):
    if x0 > 100:
        return 0.0
    res, _ = quad(lambda u: u**2 / np.expm1(u), x0, max(x0 + 80.0, 100.0),
                  limit=200)
    return res

def phi_exact(r_ld, nu_threshold=None):
    """Ionizing photon flux of the local blackbody above ``nu_threshold``
    (default: 1 Ryd). A higher threshold mimics the gravitational
    redshift cut for photons that must climb out from radius r."""
    no_torque = np.sqrt(max(0.0, (r_ld - r_in) / r_ld))
    T = TV1 * (r_ld / R0)**(-0.75) * no_torque
    if T < 100:
        return 0.0
    nu0 = NU_0 if nu_threshold is None else nu_threshold
    x0 = H_PL * nu0 / (K_B * T)
    pref = 2.0 * np.pi / C_CGS**2 * (K_B / H_PL)**3
    return pref * T**3 * wien_integral(x0)

Q_disc, _ = quad(lambda rl: phi_exact(rl) * 2.0 * np.pi * rl * LD_CM**2,
                 r_in, 1.0, limit=400)
print(f'Q_lamp = {Q_lamp:.3e}, Q_disc = {Q_disc:.3e}, '
      f'ratio = {Q_disc / Q_lamp:.3e}')

# ---- Radial profiles ----
R_PLOT_MAX = 30.0
r = np.logspace(np.log10(r_in * 1.05), np.log10(R_PLOT_MAX), 500)

# Lamp flux on a face-on element (cos(theta_inc) cancels in the ratio;
# use the flat-disc geometry for the top panel)
d = np.sqrt(r**2 + HLAMP**2)
cos_inc = HLAMP / d
Phi_lamp = Q_lamp * cos_inc / (4.0 * np.pi * (d * LD_CM)**2)

# Local viscous ionizing flux (exact Planck integral), same as Fig A1
Phi_visc_local = np.array([phi_exact(rl) for rl in r])

# ---- Emitting inner-disc grid (for the exact self-irradiation integral) ----
# Emission beyond ~1 ld is negligible (see Fig A1). Lambertian local
# intensity I = Phi_exact(r_e)/pi photons s^-1 cm^-2 sr^-1.
N_RE, N_PHI = 200, 720
r_e = np.logspace(np.log10(r_in * 1.001), 0.0, N_RE)          # ld
phi_e = np.linspace(0.0, 2.0 * np.pi, N_PHI, endpoint=False)
Phi_e = np.array([phi_exact(rl) for rl in r_e])               # cm^-2 s^-1
# Area elements r_e dr_e dphi (in cm^2)
dr_e = np.gradient(r_e)
dA = (r_e * dr_e)[:, None] * (2.0 * np.pi / N_PHI) * LD_CM**2  # (N_RE, 1)
x_e = (r_e[:, None] * np.cos(phi_e)[None, :])                  # (N_RE, N_PHI) ld
y_e = (r_e[:, None] * np.sin(phi_e)[None, :])

# ---- Exact Schwarzschild geodesic kernel (near + far side) ----
# For each emitting radius R (in GM/c^2; r_s = 2) a fan of escaping
# rays is integrated (u'' + u = 3u^2) to build the exact map
# alpha(psi): local emission angle (static observer; tan(alpha) =
# sqrt(1-2/R) r dphi/dr) vs total sweep angle psi to infinity. psi is
# the angle at the BH between (BH -> emitter) and the asymptotic
# (BH -> rim) direction; psi -> 180 deg covers far-side over-the-top
# (returning) rays. Capture is automatic: captured rays never enter
# the escaping fan. Secondary images (sweep 2pi - psi, passing under
# the disc plane) are blocked by the opaque disc and ignored.
#
# Kernel modification relative to flat space, per (R, psi):
#   M(R, psi) = [sin(alpha)/sin(psi)] * [d cos(alpha)/d cos(psi)]
# (local Lambertian factor ratio x lensing Jacobian); the bent
# integrand uses cos_e * Jac = s_t * M with s_t = z_t/rho_t.
#
# Gravitational redshift (photon-number effects):
#   - rate dilution sqrt(1 - 2/R) (time dilation of the emitter)
#   - ionization threshold blueshift: only photons EMITTED above
#     nu_0/sqrt(1-2/R) are still ionizing at infinity.
#
# Annuli inside the photon sphere (R < 3 GM/c^2) are excluded from the
# bent integral (narrow escape cones, near-total redshift); their flat
# Newtonian contribution is reported separately. NB: geodesics are
# Schwarzschild while the disc model assumes a spinning BH
# (r_in = 2 GM/c^2); treat the R < 3 zone as indicative only.
RS_LD = 2.0 * r_g_ld
COS_DPHI = np.cos(phi_e)[None, :]                              # (1, N_PHI)
R_GM = r_e / r_g_ld                                            # emitter radii in GM/c^2
PH_SPHERE = 3.0
GEO_OK = (R_GM > PH_SPHERE + 0.05)[:, None]                    # (N_RE, 1)

PSI_GRID = np.linspace(0.005, np.pi - 0.005, 720)


def build_alpha_psi_fan(R_gm, n_alpha=160):
    """Integrate a fan of escaping rays from radius R_gm (GM/c^2 units).

    Returns (psi_fan, alpha_fan) sorted in psi, for the primary image.
    """
    # Escape window: all outgoing rays (alpha < 90 deg) escape for
    # R > photon sphere; ingoing rays escape while b > sqrt(27), i.e.
    # alpha < pi - arcsin(sqrt(27(1-2/R))/R).
    sin_crit = np.sqrt(27.0 * (1.0 - 2.0 / R_gm)) / R_gm
    sin_crit = min(sin_crit, 1.0)
    a_max = np.pi - np.arcsin(sin_crit)
    # Dense fan, packed toward a_max where psi diverges (winding)
    a_lin = np.linspace(0.005, a_max - 0.05, n_alpha // 2)
    a_dense = a_max - np.geomspace(0.05, 1e-5, n_alpha // 2)
    alpha_fan = np.concatenate([a_lin, a_dense])

    u0 = 1.0 / R_gm
    u = np.full_like(alpha_fan, u0)
    w = -u0 * np.sqrt(1.0 - 2.0 / R_gm) / np.tan(alpha_fan)
    b_imp = R_gm * np.sin(alpha_fan) / np.sqrt(1.0 - 2.0 / R_gm)

    U_STOP = 1.0 / 2000.0
    DCHI = 0.005
    chi = 0.0
    psi_fan = np.full_like(alpha_fan, np.nan)
    active = np.ones(alpha_fan.shape, dtype=bool)

    def deriv(u_, w_):
        return w_, 3.0 * u_**2 - u_

    while active.any() and chi < 2.5 * np.pi:
        # RK4 step for active rays
        ua, wa = u[active], w[active]
        k1u, k1w = deriv(ua, wa)
        k2u, k2w = deriv(ua + 0.5 * DCHI * k1u, wa + 0.5 * DCHI * k1w)
        k3u, k3w = deriv(ua + 0.5 * DCHI * k2u, wa + 0.5 * DCHI * k2w)
        k4u, k4w = deriv(ua + DCHI * k3u, wa + DCHI * k3w)
        u_new = ua + DCHI / 6.0 * (k1u + 2 * k2u + 2 * k3u + k4u)
        w_new = wa + DCHI / 6.0 * (k1w + 2 * k2w + 2 * k3w + k4w)
        chi += DCHI
        # Escape: u dropped below U_STOP -> record sweep + asymptote corr.
        esc = u_new < U_STOP
        idx_active = np.where(active)[0]
        if esc.any():
            idx_esc = idx_active[esc]
            # asymptotic residual sweep ~ b*u at the stop radius
            psi_fan[idx_esc] = chi + b_imp[idx_esc] * u_new[esc]
        # Capture safeguard (should not trigger inside the window)
        cap = u_new > 1.0 / 2.02
        u[idx_active] = u_new
        w[idx_active] = w_new
        active[idx_active[esc | cap]] = False

    ok = np.isfinite(psi_fan)
    psi_f, alpha_f = psi_fan[ok], alpha_fan[ok]
    order = np.argsort(psi_f)
    return psi_f[order], alpha_f[order]


def build_kernel_table():
    """M(R, psi) on (R_KERNEL, PSI_GRID), interpolated to r_e rows."""
    R_KERNEL = np.geomspace(PH_SPHERE + 0.06, R_GM.max() * 1.01, 30)
    M_tab = np.zeros((R_KERNEL.size, PSI_GRID.size))
    for i, Rg in enumerate(R_KERNEL):
        psi_f, alpha_f = build_alpha_psi_fan(Rg)
        # primary image: restrict to psi < pi (monotone branch)
        keep = psi_f < np.pi
        psi_f, alpha_f = psi_f[keep], alpha_f[keep]
        alpha_on_grid = np.interp(PSI_GRID, psi_f, alpha_f,
                                  left=np.nan, right=np.nan)
        cos_a = np.cos(alpha_on_grid)
        jac = np.gradient(cos_a, np.cos(PSI_GRID))
        M = np.sin(alpha_on_grid) / np.sin(PSI_GRID) * np.abs(jac)
        M_tab[i] = np.where(np.isfinite(M), M, 0.0)
    # interpolate kernel rows onto the emitting grid radii
    M_full = np.zeros((r_e.size, PSI_GRID.size))
    for j in range(PSI_GRID.size):
        M_full[:, j] = np.interp(R_GM, R_KERNEL, M_tab[:, j],
                                 left=0.0, right=M_tab[-1, j])
    M_full[~GEO_OK[:, 0], :] = 0.0
    return M_full

print('Building exact geodesic kernel (fan integration)...')
M_FULL = build_kernel_table()                                  # (N_RE, N_PSI)
print('  kernel done.')

# Redshifted emissivity: rate dilution + blueshifted threshold
REDSHIFT_FAC = np.sqrt(np.maximum(1.0 - 2.0 / R_GM, 0.0))
Phi_e_red = np.array([
    g * phi_exact(rl, nu_threshold=NU_0 / g) if g > 0 else 0.0
    for rl, g in zip(r_e, REDSHIFT_FAC)
])

def phi_visc_self(r_t, z_t, slope_t, bending=False):
    """Exact self-irradiation photon flux at bowl point (r_t, z_t).

    Integrates the Lambertian emission of the flat inner disc (z=0)
    over (r_e, phi_e). Target normal n = (-slope, 0, 1)/sqrt(1+slope^2).

    With ``bending=True`` the emission cosine and solid angle are
    modified by the Beloborodov (2002) direct-ray prescription (see
    above). Returns (F_direct, F_skipped_newtonian): the second term is
    the flat-space contribution of annuli inside the validity limit
    R <= 2 r_s (zero when bending=False, where it is included in
    F_direct).
    """
    dx = x_e - r_t          # ld
    dy = y_e
    dz = -z_t
    d2_ld = dx**2 + dy**2 + dz**2
    d_ld = np.sqrt(d2_ld)
    # Flat-space emission angle from the (vertical) emitter normal
    cos_e_flat = z_t / d_ld
    # Incidence on tilted target surface:
    # n = (-slope, 0, 1)/nrm; unit vector target->emitter = (dx, dy, dz)/d
    nrm = np.sqrt(1.0 + slope_t**2)
    cos_inc = (-slope_t * dx + dz) / (d_ld * nrm)

    if not bending:
        good = (cos_e_flat > 0) & (cos_inc > 0)
        integrand = np.where(
            good,
            Phi_e[:, None] / np.pi * cos_e_flat * cos_inc
            / (d2_ld * LD_CM**2),
            0.0)
        return float(np.sum(integrand * dA)), 0.0

    # --- Exact Schwarzschild bending (near + far side, primary image) ---
    # Asymptotic (rim) direction seen from the BH; emitter at
    # ~light-minutes, rim at ~light-days, so psi is R-independent:
    #   cos(psi) = r_t cos(dphi) / sqrt(r_t^2 + z_t^2)
    rho_t = np.sqrt(r_t**2 + z_t**2)
    s_t = z_t / rho_t
    cos_psi = (r_t / rho_t) * COS_DPHI                  # (1, N_PHI)
    psi = np.arccos(np.clip(cos_psi, -1.0, 1.0))
    # Look up M(R, psi) by linear interpolation on the uniform PSI_GRID
    dpsi = PSI_GRID[1] - PSI_GRID[0]
    f_idx = (psi - PSI_GRID[0]) / dpsi
    i0 = np.clip(f_idx.astype(int), 0, PSI_GRID.size - 2)
    wgt = np.clip(f_idx - i0, 0.0, 1.0)
    i0 = np.broadcast_to(i0, (r_e.size,) + i0.shape[1:])
    wgt = np.broadcast_to(wgt, i0.shape)
    rows = np.arange(r_e.size)[:, None]
    M = (1.0 - wgt) * M_FULL[rows, i0] + wgt * M_FULL[rows, i0 + 1]
    # Bent kernel: cos_e * Jacobian = s_t * M; redshifted emissivity
    good = GEO_OK & (M > 0) & (cos_inc > 0)
    integrand = np.where(
        good,
        Phi_e_red[:, None] / np.pi * (s_t * M) * cos_inc
        / (d2_ld * LD_CM**2),
        0.0)
    F_bend = float(np.sum(integrand * dA))

    # Flat Newtonian contribution of annuli inside the photon sphere
    # (excluded from the geodesic kernel), reported separately.
    good_skip = (~GEO_OK) & (cos_e_flat > 0) & (cos_inc > 0)
    integrand_skip = np.where(
        good_skip,
        Phi_e[:, None] / np.pi * cos_e_flat * cos_inc
        / (d2_ld * LD_CM**2),
        0.0)
    F_skip = float(np.sum(integrand_skip * dA))
    return F_bend, F_skip

# ---- Per-h1 surface fluxes ----
curves = {}
curves_bend = {}
for h1 in H1_LIST:
    z = h1 * (r / R0)**BETA
    slope = BETA * z / r            # dz/dr
    norm = np.sqrt(1.0 + slope**2)

    # Lamp on the tilted bowl surface:
    # normal n = (-slope, 1)/norm; direction to lamp l = (-r, hlamp - z)/d
    d_lamp = np.sqrt(r**2 + (HLAMP - z)**2)
    cos_lamp = (slope * r + (HLAMP - z)) / (d_lamp * norm)
    cos_lamp = np.clip(cos_lamp, 0.0, 1.0)
    Phi_lamp_surf = Q_lamp * cos_lamp / (4.0 * np.pi * (d_lamp * LD_CM)**2)

    # Exact inner-disc self-irradiation integral at every target radius
    F_disc = np.array([
        phi_visc_self(r_t, z_t, s_t)[0] for r_t, z_t, s_t in zip(r, z, slope)
    ])
    res_bend = [phi_visc_self(r_t, z_t, s_t, bending=True)
                for r_t, z_t, s_t in zip(r, z, slope)]
    F_disc_bend = np.array([fr[0] for fr in res_bend])
    F_skip_newton = np.array([fr[1] for fr in res_bend])

    # No disc surface beyond r_out: truncate all surface fluxes there
    beyond = r > R_OUT
    Phi_lamp_surf = np.where(beyond, np.nan, Phi_lamp_surf)
    F_disc = np.where(beyond, np.nan, F_disc)
    F_disc_bend = np.where(beyond, np.nan, F_disc_bend)
    F_skip_newton = np.where(beyond, np.nan, F_skip_newton)

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(Phi_lamp_surf > 0, F_disc / Phi_lamp_surf, np.nan)
        ratio_bend = np.where(Phi_lamp_surf > 0, F_disc_bend / Phi_lamp_surf,
                              np.nan)
    curves[h1] = (Phi_lamp_surf, F_disc, ratio)
    curves_bend[h1] = (Phi_lamp_surf, F_disc_bend, ratio_bend)
    imax = np.nanargmax(ratio)
    print(f'H(r_out)={h1}: max F_disc/F_lamp = {np.nanmax(ratio):.3e} '
          f'(exact GR near+far: {np.nanmax(ratio_bend):.3e}) '
          f'at r = {r[imax]:.2f} ld; '
          f'photon-sphere annuli (R<3) flat share there = '
          f'{F_skip_newton[imax] / F_disc[imax]:.2f} of flat F_disc')


def make_figure(outname, xmin, xmax, xscale, ylim1, ylim2,
                show_rmin=True, legend1_loc='lower left',
                curve_set=None, self_label=r'$\Phi_{\rm visc,self}$'):
    if curve_set is None:
        curve_set = curves
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True,
                                   gridspec_kw={'hspace': 0.05})

    ax1.plot(r, np.log10(np.maximum(Phi_visc_local, 1e-30)), color='black',
             lw=2.5, ls='-', alpha=0.85,
             label=r'$\Phi_{\rm visc}~(\rm local,~exact)$')

    for h1, col in zip(H1_LIST, COLORS):
        Phi_lamp_surf, F_disc, ratio = curve_set[h1]
        lab = rf'$H(r_{{\rm out}}) = {h1}~\rm ld$'
        with np.errstate(invalid='ignore'):
            ax1.plot(r, np.log10(np.maximum(Phi_lamp_surf, 1e-30)), color=col,
                     lw=2, ls='-', alpha=0.9)
            ax1.plot(r, np.log10(np.maximum(F_disc, 1e-30)), color=col, lw=3,
                     ls='--', alpha=0.9, label=lab)  # Phi_visc,self
            F_total = Phi_lamp_surf + Phi_visc_local + F_disc
            ax1.plot(r, np.log10(F_total), color=col, lw=2, ls=':', alpha=0.9)
        ax2.plot(r, ratio, color=col, lw=3, ls='--', alpha=0.9, label=lab)

    # line-style key: solid = lamp on bowl, dashed = disc irrad, dotted = total
    ax1.plot([], [], color='gray', lw=2, ls='-',
             label=r'$\Phi_{\rm irrad}~(\rm lamp~on~bowl)$')
    ax1.plot([], [], color='gray', lw=3, ls='--', label=self_label)
    ax1.plot([], [], color='gray', lw=2, ls=':',
             label=r'$\Phi_{\rm total}$')

    # Panel 1 styling
    ax1.set_ylim(*ylim1)
    ax1.set_ylabel(r'$\log\Phi~\rm [photons~s^{-1}~cm^{-2}]$')
    if show_rmin:
        ax1.axvline(R_PILLAR_MIN, color='purple', ls='-.', lw=2, alpha=0.8)
        ax1.text(R_PILLAR_MIN * 0.9, ylim1[0] + 2.0,
                 r'$r_{\rm min}^{\rm pillar}$',
                 color='purple', fontsize=15, rotation=90,
                 ha='right', va='bottom')
    ax1.axvline(R_OUT, color='purple', ls='-.', lw=2, alpha=0.8)
    if xscale == 'log':
        x_rout_text = R_OUT * 1.07
    else:
        x_rout_text = R_OUT + 0.015 * (xmax - xmin)
    ax1.text(x_rout_text, ylim1[0] + 2.0, r'$r_{\rm out}$',
             color='purple', fontsize=15, rotation=90,
             ha='left', va='bottom')
    ax1.legend(fontsize=14, loc=legend1_loc, framealpha=0.95)

    # Panel 2 styling
    ax2.axhline(1.0, color='gray', ls=':', lw=2, alpha=0.8)
    ax2.set_yscale('log')
    ax2.set_ylim(*ylim2)
    ax2.set_xscale(xscale)
    ax2.set_xlim(xmin, xmax)
    ax2.set_xlabel(r'$\rm Radius~[ld]$')
    ax2.set_ylabel(self_label + r'$\,/\,\Phi_{\rm lamp}$')
    if show_rmin:
        ax2.axvline(R_PILLAR_MIN, color='purple', ls='-.', lw=2, alpha=0.8)
    ax2.axvline(R_OUT, color='purple', ls='-.', lw=2, alpha=0.8)
    ax2.legend(fontsize=14, loc='upper left', framealpha=0.95)

    for ax in (ax1, ax2):
        ax.tick_params(which='major', length=8, width=1.5, direction='in',
                       top=True, right=True)
        ax.tick_params(which='minor', length=4, width=1, direction='in',
                       top=True, right=True)
        ax.minorticks_on()

    os.makedirs('plots', exist_ok=True)
    plt.savefig(outname, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {outname}')


# Full view (log x)
make_figure('plots/fig_self_irradiation.png',
            r_in, R_PLOT_MAX, 'log', (8, 26), (1e-14, 3))

# Zoom on the rim (linear x)
make_figure('plots/fig_self_irradiation_rim.png',
            16.0, 20.5, 'linear', (14, 22), (1e-6, 3), show_rmin=False,
            legend1_loc='upper left')

# Light bending: exact Schwarzschild geodesics (near + far side)
BEND_LABEL = r'$\Phi_{\rm visc,self}^{\rm GR}$'
make_figure('plots/fig_self_irradiation_bend.png',
            r_in, R_PLOT_MAX, 'log', (8, 26), (1e-14, 3),
            curve_set=curves_bend, self_label=BEND_LABEL)
make_figure('plots/fig_self_irradiation_bend_rim.png',
            16.0, 20.5, 'linear', (14, 22), (1e-6, 3), show_rmin=False,
            legend1_loc='upper left',
            curve_set=curves_bend, self_label=BEND_LABEL)
