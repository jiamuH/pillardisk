r"""Back-trace rays from the rim (r=20 ld, z=1 ld) to the inner disc.

Rim direction seen from the BH: elevation eta = arctan(1/20) = 2.862 deg,
at distance ~5012 GM/c^2 (1 ld = 250.6 GM/c^2 for M_BH = 7e7 Msun).

For each emitting radius R we shoot Schwarzschild geodesics
(u'' + u = 3u^2, GM/c^2 units) and bisect on the local emission angle
alpha (from outward radial) so the asymptotic direction matches the rim:
  - near side (dphi = 0):  asymptotic theta = eta
  - far side  (dphi = pi): exit theta = eta after bending over the BH

Also prints the per-emitter kernel modification vs flat space:
  intensity ratio sin(alpha_bent)/sin(alpha_flat),
  Jacobian d cos(alpha)/d cos(psi) (numerical),
  net = product; compared to Beloborodov (1 - r_s/R)^{3/2}.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
    'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

R_HOR = 2.0
R_PH = 3.0
R_FAR = 3000.0           # "infinity" for the shooting map

LD_GM = 250.6            # 1 ld in GM/c^2 for M_BH = 7e7 Msun
RIM_R = 20.0 * LD_GM
RIM_Z = 1.0 * LD_GM
ETA = np.degrees(np.arctan2(RIM_Z, RIM_R))   # 2.862 deg


def trace(r0, alpha_deg, theta0=0.0, sweep_sign=+1.0, chi_max=4 * np.pi,
          r_stop=R_FAR):
    a = np.radians(alpha_deg)
    u0 = 1.0 / r0
    # Local angle of a static observer: tan(alpha) = sqrt(1-rs/r) r dphi/dr
    # (proper radial length is dr/sqrt(1-rs/r)), so
    #   du/dchi = -u sqrt(1-rs/r) / tan(alpha).
    w0 = -u0 * np.sqrt(1.0 - R_HOR / r0) / np.tan(a)

    def rhs(chi, y):
        return [y[1], 3.0 * y[0]**2 - y[0]]

    def horizon(chi, y):
        return 1.0 / max(y[0], 1e-14) - R_HOR * 1.02
    horizon.terminal = True
    horizon.direction = -1

    def far(chi, y):
        return 1.0 / max(y[0], 1e-14) - r_stop
    far.terminal = True
    far.direction = 1

    sol = solve_ivp(rhs, (0, chi_max), [u0, w0], rtol=1e-11, atol=1e-13,
                    events=[horizon, far], max_step=0.02)
    r_arr = 1.0 / np.maximum(sol.y[0], 1e-14)
    theta = theta0 + sweep_sign * sol.t
    captured = bool(sol.t_events[0].size)
    return r_arr, theta, captured


def exit_angle(r0, alpha_deg, theta0, sweep_sign):
    """Polar angle (deg, in (-180, 180]) of the escape direction."""
    r_arr, theta, captured = trace(r0, alpha_deg, theta0, sweep_sign)
    if captured:
        return np.nan
    x = r_arr[-1] * np.cos(theta[-1])
    y = r_arr[-1] * np.sin(theta[-1])
    return np.degrees(np.arctan2(y, x))


def bisect(f, lo, hi, tol=1e-9):
    flo, fhi = f(lo), f(hi)
    assert np.isfinite(flo) and np.isfinite(fhi) and flo * fhi < 0, \
        f'no bracket: f({lo})={flo}, f({hi})={fhi}'
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if not np.isfinite(fm):
            # captured: move toward the escaping side (assume lo escapes)
            hi = mid
            continue
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


# ---- Near side: solve theta_inf(alpha) = ETA for several R ----
print(f'rim elevation eta = {ETA:.3f} deg')
near = {}
for R0 in [3.3, 6.0, 10.0]:
    a_sol = bisect(lambda a: exit_angle(R0, a, 0.0, +1.0) - ETA, 0.2, 30.0)
    near[R0] = a_sol
    # numerical Jacobian d cos(alpha)/d cos(psi)
    da = 0.05
    psi1 = exit_angle(R0, a_sol - da, 0.0, +1.0)
    psi2 = exit_angle(R0, a_sol + da, 0.0, +1.0)
    dcos_a = np.cos(np.radians(a_sol - da)) - np.cos(np.radians(a_sol + da))
    dcos_p = np.cos(np.radians(psi1)) - np.cos(np.radians(psi2))
    jac = dcos_a / dcos_p
    a_flat = ETA
    ratio_I = np.sin(np.radians(a_sol)) / np.sin(np.radians(a_flat))
    net = ratio_I * jac
    belo = (1.0 - 2.0 / R0)**1.5
    print(f'  near R={R0:5.1f}: alpha_bent={a_sol:6.3f} deg '
          f'(flat {a_flat:.3f}); I-ratio={ratio_I:.3f}, J={jac:.3f}, '
          f'net={net:.3f}  [Beloborodov (1-rs/R)^1.5 = {belo:.3f}]')

# ---- Far side: solve exit(alpha) = ETA (over the top) for R = 3.3 ----
a_far = bisect(lambda a: exit_angle(3.3, a, np.pi, -1.0) - ETA, 84.0, 86.0)
print(f'  far  R=  3.3: alpha_bent={a_far:.3f} deg from radial '
      f'(local cos(theta_e) = {np.sin(np.radians(a_far)):.3f}: '
      f'NOT grazing-suppressed)')

# ---- Figure: zoom near the BH with the back-traced rays ----
fig, ax = plt.subplots(figsize=(13, 8))
th = np.linspace(0, 2 * np.pi, 200)
ax.fill(R_HOR * np.cos(th), R_HOR * np.sin(th), color='black', zorder=5)
ax.plot(R_PH * np.cos(th), R_PH * np.sin(th), color='gray', ls=':',
        lw=1.8, zorder=4)
for sgn in (+1, -1):
    ax.plot([sgn * R_HOR, sgn * 14.5], [0, 0], color='saddlebrown', lw=5,
            solid_capstyle='butt', alpha=0.9, zorder=3)

cols = {3.3: 'royalblue', 6.0: 'mediumseagreen', 10.0: 'darkviolet'}
for R0, a_sol in near.items():
    r_arr, theta, _ = trace(R0, a_sol, 0.0, +1.0, r_stop=80.0)
    ax.plot(r_arr * np.cos(theta), r_arr * np.sin(theta), color=cols[R0],
            lw=2.5, alpha=0.9, zorder=6,
            label=rf'$\rm near,~R={R0}~GM/c^2$')
    ax.plot(R0, 0, marker='o', color=cols[R0], ms=11, zorder=7)

r_arr, theta, _ = trace(3.3, a_far, np.pi, -1.0, r_stop=80.0)
ax.plot(r_arr * np.cos(theta), r_arr * np.sin(theta), color='orangered',
        lw=2.5, alpha=0.9, zorder=6,
        label=r'$\rm far,~R=3.3~GM/c^2~(returning)$')
ax.plot(-3.3, 0, marker='o', color='orangered', ms=11, zorder=7)

# Rim direction guide
tt = np.linspace(0, 30, 5)
ax.plot(tt, tt * np.tan(np.radians(ETA)), color='black', ls='--', lw=1.5,
        alpha=0.6, zorder=2)
ax.annotate(r'$\rm to~rim:~r=20~ld,~z=1~ld~(\eta = 2.86^{\circ})$',
            xy=(24.0, 1.35), xytext=(10.0, 4.2), fontsize=16,
            arrowprops=dict(arrowstyle='->', color='black', lw=1.6))

ax.text(0.0, 3.3, r'$\rm photon~sphere$', color='gray', fontsize=14,
        ha='center', va='bottom')
ax.text(0.0, -0.95, r'$\rm BH$', color='white', fontsize=15, ha='center',
        zorder=8)

ax.set_xlim(-12, 26)
ax.set_ylim(-3.0, 8.5)
ax.set_aspect('equal')
ax.set_xlabel(r'$x~[GM/c^2]$')
ax.set_ylabel(r'$z~[GM/c^2]$')
ax.legend(fontsize=14, loc='upper left', framealpha=0.9)
ax.tick_params(which='major', length=8, width=1.5, direction='in',
               top=True, right=True)
ax.tick_params(which='minor', length=4, width=1, direction='in',
               top=True, right=True)
ax.minorticks_on()

os.makedirs('plots', exist_ok=True)
out = 'plots/fig_backtrace_rays.png'
plt.savefig(out, dpi=200, bbox_inches='tight')
print(f'wrote {out}')
