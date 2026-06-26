"""Plot Phi_irrad (lamp) and Phi_visc (Wien-tail) vs radius.

Reproduces fig_ionflux_profile.png in the paper draft. The viscous Wien-tail
expression here is the leading-order form
  Phi_visc ≈ (2*pi*nu_0^2/c^2) * (k_B T / h) * exp(-h nu_0 / (k_B T))
matching pillar_disk._compute_direct_ionizing_flux.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
    'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

# ---- Physical constants (cgs) ----
H_PL = 6.626e-27
K_B = 1.381e-16
C_CGS = 2.998e10
NU_0 = 13.6 * 1.602e-12 / H_PL  # 1 Ryd frequency (Hz)
G_CGS = 6.674e-8
M_SUN = 1.989e33
LD_CGS = 2.998e10 * 86400.0

# ---- Parameters (match config_line.yaml) ----
M_BH = 7.0e7
R_ISCO_RG = 2.0
HLAMP = 0.05            # ld
R_OUT = 20.0            # ld (disc outer radius, config_line.yaml)
ROUT = 120.0            # ld (plot extent)
LOG_PHI_INNER = 25.0    # at r_in on the flat disc (used for Q calibration)
R_PILLAR_MIN = 2.0      # ld — inner edge of pillar zone (config_line.yaml rmin)

# Temperature profile T_v(r) = tv1 (r/r0)^(-3/4)
TV1 = 446.0
R0 = 20.0
ALPHA_T = 0.75

r_g_cm = G_CGS * M_BH * M_SUN / C_CGS**2
r_g_ld = r_g_cm / LD_CGS
r_in = R_ISCO_RG * r_g_ld

print(f'r_g = {r_g_ld:.5f} ld; r_in = {r_in:.5f} ld')

# Radial grid
r = np.logspace(np.log10(r_in), np.log10(ROUT), 600)  # ld

# ---- Lamp ionizing flux (face-on, flat disc, phi=0) ----
d_in = np.sqrt(r_in**2 + HLAMP**2)
cos_theta_in = HLAMP / d_in
Q = 10.0**LOG_PHI_INNER * 4.0 * np.pi * d_in**2 / cos_theta_in

d = np.sqrt(r**2 + HLAMP**2)
cos_theta = HLAMP / d
Phi_lamp = Q * cos_theta / (4.0 * np.pi * d**2)

# ---- Viscous Wien-tail ionizing flux ----
# T_visc with Shakura-Sunyaev no-torque factor sqrt((r - rin)/r)
no_torque = np.sqrt(np.maximum(0.0, (r - r_in) / r))
T_visc = TV1 * (r / R0)**(-ALPHA_T) * no_torque
A = 2.0 * np.pi * NU_0**2 / C_CGS**2 * (K_B / H_PL)
x = H_PL * NU_0 / (K_B * np.maximum(T_visc, 100.0))
x = np.clip(x, 0.0, 500.0)
Phi_visc = A * T_visc * np.exp(-x)

# KDH's alternative form (R-J + Wien hybrid; see notes):
#   Phi = (4*pi/c^2) (kT/h)^3 (x0 + 1) exp(-x0)
A_kdh = 4.0 * np.pi / C_CGS**2 * (K_B / H_PL)**3
Phi_visc_kdh = A_kdh * T_visc**3 * (x + 1.0) * np.exp(-x)

# Exact numerical integral
#   Phi_exact = (2*pi/c^2) (kT/h)^3 \int_{x0}^infty x^2/(e^x - 1) dx
def _wien_integral(x0):
    if x0 > 100:
        return 0.0
    res, _ = quad(lambda u: u**2 / np.expm1(u), x0, max(x0 + 80.0, 100.0),
                  limit=200)
    return res

prefac = 2.0 * np.pi / C_CGS**2 * (K_B / H_PL)**3
Phi_visc_exact = np.array([
    prefac * T_v**3 * _wien_integral(H_PL * NU_0 / (K_B * max(T_v, 100.0)))
    for T_v in T_visc
])

Phi_total = Phi_lamp + Phi_visc_exact

i_tmax = np.argmax(T_visc)
T_visc_max = T_visc[i_tmax]
print(f'T_visc,max = {T_visc_max:.0f} K at r = {r[i_tmax]:.4f} ld')

# ---- Plot ----
# semilog x; y is log10(Phi)
fig, ax = plt.subplots(figsize=(10, 8))

logL = np.log10(np.maximum(Phi_lamp, 1e-30))
logV = np.log10(np.maximum(Phi_visc, 1e-30))
logT = np.log10(np.maximum(Phi_total, 1e-30))

ax.semilogx(r, logL, color='royalblue', lw=3, alpha=0.9,
            label=r'$\Phi_{\rm irrad}~(\rm lamp)$')
ax.semilogx(r, logV, color='orangered', lw=3, ls='--', alpha=0.9,
            label=r'$\Phi_{\rm visc}~(\nu_0^2\,kT/h)$')
logV_kdh = np.log10(np.maximum(Phi_visc_kdh, 1e-30))
ax.semilogx(r, logV_kdh, color='darkviolet', lw=3, ls=':', alpha=0.9,
            label=r'$\Phi_{\rm visc}~(\nu_0\,(kT/h)^2)$')
logV_exact = np.log10(np.maximum(Phi_visc_exact, 1e-30))
ax.semilogx(r, logV_exact, color='black', lw=2.5, ls='-', alpha=0.85,
            label=r'$\Phi_{\rm visc}~(\rm exact)$')
ax.semilogx(r, logT, color='black', lw=2.5, ls=':', alpha=0.9,
            label=r'$\Phi_{\rm total}$')

YMIN, YMAX = 12.0, 26.0
XMIN = 1.5e-3

# Cloudy grid band log Phi in [15, 21]
ax.axhspan(15, 21, color='forestgreen', alpha=0.12)
ax.text(0.015, 18.0, r'$\rm Cloudy~grid~range$',
        color='forestgreen', fontsize=16, ha='left', va='center',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none',
                  boxstyle='round,pad=0.3'))

# r < r_ISCO band
ax.axvspan(XMIN, r_in, color='gray', alpha=0.25)
ax.text(np.sqrt(XMIN * r_in), 0.5 * (YMIN + YMAX),
        r'$\rm r < r_{ISCO}$',
        color='dimgray', fontsize=15, rotation=90, ha='center', va='center')

# T_visc,max annotation: point to peak of Phi_visc, place text in empty upper-right
i_peak = np.argmax(Phi_visc)
T_at_peak = T_visc[i_peak]
ax.annotate(rf'$T_{{\rm visc,max}} \approx {T_at_peak:,.0f}~\rm K$',
            xy=(r[i_peak], logV[i_peak]),
            xytext=(np.sqrt(0.01 * 0.1), 23.0),
            color='orangered', fontsize=16, ha='center',
            arrowprops=dict(arrowstyle='->', color='orangered', lw=1.5),
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none',
                      boxstyle='round,pad=0.3'))

# Pillar inner edge r_min
ax.axvline(R_PILLAR_MIN, color='purple', ls='-.', lw=2, alpha=0.8)
ax.text(R_PILLAR_MIN * 0.9, 22.5,
        r'$r_{\rm min}^{\rm pillar}$',
        color='purple', fontsize=15, rotation=90, ha='right', va='center')

# Disc outer edge r_out
ax.axvline(R_OUT, color='purple', ls='-.', lw=2, alpha=0.8)
ax.text(R_OUT * 0.9, 22.5,
        r'$r_{\rm out}$',
        color='purple', fontsize=15, rotation=90, ha='right', va='center')

ax.set_xlim(XMIN, ROUT)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel(r'$\rm Radius~[ld]$')
ax.set_ylabel(r'$\log\Phi~\rm [photons~s^{-1}~cm^{-2}]$')

ax.tick_params(which='major', length=8, width=1.5, direction='in',
               top=True, right=True)
ax.tick_params(which='minor', length=4, width=1, direction='in',
               top=True, right=True)
ax.minorticks_on()

# Top axis: r_g
ax2 = ax.twiny()
ax2.set_xscale('log')
ax2.set_xlim(XMIN / r_g_ld, ROUT / r_g_ld)
ax2.set_xlabel(r'$\rm Radius~[r_g]$')
ax2.tick_params(which='major', length=8, width=1.5, direction='in')
ax2.tick_params(which='minor', length=4, width=1, direction='in')
ax2.minorticks_on()

handles, labels = ax.get_legend_handles_labels()
order = [labels.index(r'$\Phi_{\rm total}$'),
         labels.index(r'$\Phi_{\rm irrad}~(\rm lamp)$'),
         labels.index(r'$\Phi_{\rm visc}~(\rm exact)$'),
         labels.index(r'$\Phi_{\rm visc}~(\nu_0^2\,kT/h)$'),
         labels.index(r'$\Phi_{\rm visc}~(\nu_0\,(kT/h)^2)$')]
ax.legend([handles[i] for i in order], [labels[i] for i in order],
          loc='lower left', fontsize=13, framealpha=0.95)

os.makedirs('plots', exist_ok=True)
out = 'plots/fig_ionflux_profile.png'
plt.savefig(out, dpi=200, bbox_inches='tight')
print(f'wrote {out}')
