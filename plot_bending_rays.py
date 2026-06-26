r"""Schwarzschild photon trajectories illustrating the two bending cases.

Null geodesics in the orbital (figure) plane obey
    d^2 u / d chi^2 + u = 3 u^2     (units GM/c^2 = 1, so r_s = 2),
with u = 1/r and chi the sweep angle from emission. Emission at radius
r0 and local angle alpha from the OUTWARD radial direction gives
    u(0) = 1/r0,   du/dchi(0) = -u(0)/tan(alpha).

Case 1 (direct, near side): emitter at theta=0 (x=+r0), mildly bent
rays that escape to the rim side (+x). Compared against straight rays
with the same local emission angle: bending increases the sweep, so the
bent ray arrives at HIGHER asymptotic elevation; to hit the near-plane
rim the photon must be emitted more grazing (the (1-r_s/R)^{3/2}
suppression).

Case 2 (strongly bent, far side): emitter at theta=pi (x=-r0), rays
emitted up/inward that pass near the photon sphere (r=3), get bent over
the top of the BH, and exit toward +x at low elevation: the
"elevated-lamp" returning radiation. One slightly more inward ray is
captured.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
    'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

R_HOR = 2.0      # horizon (r_s) in GM/c^2
R_PH = 3.0       # photon sphere
R_MAX = 60.0


def trace_ray(r0, alpha_deg, theta0, sweep_sign, chi_max=4.0 * np.pi):
    """Integrate a null geodesic; return (x, y, fate, theta_asym).

    theta(chi) = theta0 + sweep_sign * chi. fate: 'escape'|'capture'.
    """
    a = np.radians(alpha_deg)
    u0 = 1.0 / r0
    # Local angle of a static observer: tan(alpha) = sqrt(1-rs/r) r dphi/dr
    w0 = -u0 * np.sqrt(1.0 - R_HOR / r0) / np.tan(a)

    def rhs(chi, y):
        u, w = y
        return [w, 3.0 * u**2 - u]

    def hit_horizon(chi, y):
        return 1.0 / max(y[0], 1e-12) - (R_HOR * 1.02)
    hit_horizon.terminal = True
    hit_horizon.direction = -1

    def escaped(chi, y):
        return 1.0 / max(y[0], 1e-12) - R_MAX
    escaped.terminal = True
    escaped.direction = 1

    sol = solve_ivp(rhs, (0.0, chi_max), [u0, w0], rtol=1e-10, atol=1e-12,
                    dense_output=True, events=[hit_horizon, escaped],
                    max_step=0.02)
    chi = sol.t
    u = sol.y[0]
    r_arr = 1.0 / np.maximum(u, 1e-12)
    theta = theta0 + sweep_sign * chi
    x = r_arr * np.cos(theta)
    y = r_arr * np.sin(theta)
    fate = 'capture' if sol.t_events[0].size else 'escape'
    return x, y, fate, np.degrees(theta[-1] % (2 * np.pi))


fig, ax = plt.subplots(figsize=(13, 8))

# --- Black hole, photon sphere, disc plane ---
th = np.linspace(0, 2 * np.pi, 200)
ax.fill(R_HOR * np.cos(th), R_HOR * np.sin(th), color='black', zorder=5)
ax.plot(R_PH * np.cos(th), R_PH * np.sin(th), color='gray', ls=':',
        lw=1.8, zorder=4)
# Disc plane (emitting inner disc, peak annuli near r ~ 3.3)
for sgn in (+1, -1):
    ax.plot([sgn * R_HOR, sgn * 14.5], [0, 0], color='saddlebrown', lw=5,
            solid_capstyle='butt', alpha=0.9, zorder=3)

# --- Case 1: direct rays, near-side emitter at x = +6 ---
R0_NEAR = 6.0
col1 = 'royalblue'
for alpha in [8.0, 20.0, 35.0]:
    x, y, fate, th_asym = trace_ray(R0_NEAR, alpha, 0.0, +1.0)
    ax.plot(x, y, color=col1, lw=2.5, alpha=0.9, zorder=6)
    print(f'near-side alpha={alpha:5.1f}: fate={fate}, '
          f'asymptotic theta={th_asym:.1f} deg')
# Straight (flat-space) comparison for the middle ray
a = np.radians(20.0)
t = np.linspace(0, 40, 10)
ax.plot(R0_NEAR + t * np.cos(a), t * np.sin(a), color=col1, lw=1.8,
        ls='--', alpha=0.7, zorder=6)
ax.plot(R0_NEAR, 0, marker='o', color=col1, ms=12, zorder=7)

# --- Case 2: strongly bent rays, far-side emitter at x = -3.3 ---
# Escape over the BH requires impact parameter b = r0 sin(a)/sqrt(1-2/r0)
# > sqrt(27); for r0 = 3.3 that is alpha in (90, ~99.2) deg.
R0_FAR = 3.3
col2 = 'orangered'
for alpha in [82.0, 84.5]:
    x, y, fate, th_asym = trace_ray(R0_FAR, alpha, np.pi, -1.0)
    ax.plot(x, y, color=col2, lw=2.5, alpha=0.9, zorder=6)
    print(f'far-side  alpha={alpha:5.1f}: fate={fate}, '
          f'asymptotic theta={th_asym:.1f} deg')
# A captured ray (slightly more inward)
x, y, fate, th_asym = trace_ray(R0_FAR, 100.0, np.pi, -1.0)
ax.plot(x, y, color='gray', lw=2.0, ls='-', alpha=0.8, zorder=6)
print(f'far-side  alpha=100.0: fate={fate}')
ax.plot(-R0_FAR, 0, marker='o', color=col2, ms=12, zorder=7)

# --- Annotations ---
ax.annotate(r'$\rm to~rim~(\sim 5000~GM/c^2)$',
            xy=(14.0, 1.8), xytext=(7.0, 4.6),
            fontsize=17, color='black',
            arrowprops=dict(arrowstyle='->', color='black', lw=1.8))
ax.text(6.2, -1.1, r'$\rm direct~(near~side)$', color=col1, fontsize=17,
        ha='left')
ax.text(-11.6, 2.6, r'$\rm strongly~bent$', color=col2, fontsize=17,
        ha='left')
ax.text(-11.6, 1.7, r'$\rm (far~side,~returning)$', color=col2,
        fontsize=17, ha='left')
ax.text(-2.1, -2.6, r'$\rm captured$', color='gray', fontsize=15,
        ha='center')
ax.text(0.0, 3.35, r'$\rm photon~sphere~(3\,GM/c^2)$', color='gray',
        fontsize=14, ha='center', va='bottom')
ax.text(0.0, -0.95, r'$\rm BH$', color='white', fontsize=15, ha='center',
        zorder=8)
ax.text(11.0, -0.85, r'$\rm disc$', color='saddlebrown', fontsize=15,
        ha='center')

# Dashed straight-ray label
ax.text(10.6, 8.0, r'$\rm flat~space~(same~local~angle)$', color=col1,
        fontsize=14, ha='center', alpha=0.85)

ax.set_xlim(-12, 15)
ax.set_ylim(-3.2, 9.5)
ax.set_aspect('equal')
ax.set_xlabel(r'$x~[GM/c^2]$')
ax.set_ylabel(r'$z~[GM/c^2]$')
ax.tick_params(which='major', length=8, width=1.5, direction='in',
               top=True, right=True)
ax.tick_params(which='minor', length=4, width=1, direction='in',
               top=True, right=True)
ax.minorticks_on()

os.makedirs('plots', exist_ok=True)
out = 'plots/fig_bending_rays.png'
plt.savefig(out, dpi=200, bbox_inches='tight')
print(f'wrote {out}')
