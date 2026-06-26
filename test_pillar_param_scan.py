"""
Parameter scan: how pillar geometry affects the azimuthal-average T(r) enhancement.

We scan over h_pillar, sigma_r, sigma_phi, and h_lamp to find the combination
that maximizes the temperature boost at the pillar radius while minimizing
shadow suppression at larger radii.
"""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin

config = load_config('config_line.yaml')


def build_disk(hlamp=None):
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    if hlamp is not None:
        lamp_params['hlamp'] = hlamp
    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                    'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc',
                    'cosi', 'inclination', 'redshift', 'M_BH', 'r_isco_rg', 'f_trans',
                    'f_turb_line']
    int_params = ['nr', 'nphi']
    for d in [disk_params, temp_params, lamp_params, obs_params]:
        for k, v in d.items():
            if isinstance(v, str) and v.lower() == 'auto':
                continue
            if k in int_params:
                d[k] = int(float(v))
            elif k in float_params and v is not None:
                d[k] = float(v)
    resolve_rin(disk_params)
    if 'inclination' in obs_params:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']
    return PillarDisk(**{**disk_params, **temp_params, **lamp_params, **obs_params})


def get_T_ratio(hlamp, h_pillar, sigma_r, sigma_phi, r_p=10.0):
    """
    Compute T(r) with one pillar / T(r) without pillar.
    Returns (r, T_ratio, T_with, T_without).
    """
    disk = build_disk(hlamp=hlamp)
    # Add one pillar at phi=0
    disk.add_pillar(r_pillar=r_p, phi_pillar=0.0, height=h_pillar,
                    sigma_r=sigma_r, sigma_phi=sigma_phi,
                    modify_height=True, modify_temp=True)
    r, T_with, T_without, _, _, _ = disk.get_temperature_profile()
    T_ratio = T_with / np.maximum(T_without, 1.0)
    return r, T_ratio, T_with, T_without


def style_ax(ax):
    ax.tick_params(direction='in', which='major', length=8, width=1.8, labelsize=13,
                   top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


r_p = 10.0  # pillar radius

# ======================================================
# Scan 1: varying h_pillar/sigma_r (aspect ratio = steepness)
# ======================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Top left: vary h_pillar with fixed sigma_r
ax = axes[0, 0]
sigma_r_fix = 0.5
sigma_phi_fix = 0.2
hlamp_fix = 0.5
for h_p in [0.05, 0.1, 0.2, 0.5, 1.0]:
    r, T_ratio, _, _ = get_T_ratio(hlamp_fix, h_p, sigma_r_fix, sigma_phi_fix, r_p)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8,
                label=rf'$h_p = {h_p}~\rm ld$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.axvline(r_p, color='gray', ls=':', lw=1, alpha=0.5)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm flat}$', fontsize=16)
ax.set_title(rf'$\rm Vary~h_p~(\sigma_r={sigma_r_fix},~\sigma_\phi={sigma_phi_fix},~h_{{LP}}={hlamp_fix})$',
             fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

# Top right: vary sigma_r with fixed h_pillar
ax = axes[0, 1]
h_p_fix = 0.3
for sr in [0.1, 0.2, 0.5, 1.0, 2.0]:
    r, T_ratio, _, _ = get_T_ratio(hlamp_fix, h_p_fix, sr, sigma_phi_fix, r_p)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8,
                label=rf'$\sigma_r = {sr}~\rm ld$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.axvline(r_p, color='gray', ls=':', lw=1, alpha=0.5)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm flat}$', fontsize=16)
ax.set_title(rf'$\rm Vary~\sigma_r~(h_p={h_p_fix},~\sigma_\phi={sigma_phi_fix},~h_{{LP}}={hlamp_fix})$',
             fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

# Bottom left: vary sigma_phi
ax = axes[1, 0]
for sp in [0.05, 0.1, 0.2, 0.5, 1.0]:
    r, T_ratio, _, _ = get_T_ratio(hlamp_fix, h_p_fix, sigma_r_fix, sp, r_p)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8,
                label=rf'$\sigma_\phi = {sp}~\rm rad$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.axvline(r_p, color='gray', ls=':', lw=1, alpha=0.5)
ax.set_xlabel(r'$r~\rm [ld]$', fontsize=16)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm flat}$', fontsize=16)
ax.set_title(rf'$\rm Vary~\sigma_\phi~(h_p={h_p_fix},~\sigma_r={sigma_r_fix},~h_{{LP}}={hlamp_fix})$',
             fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

# Bottom right: vary h_lamp
ax = axes[1, 1]
for hl in [0.05, 0.1, 0.2, 0.5, 1.0]:
    r, T_ratio, _, _ = get_T_ratio(hl, h_p_fix, sigma_r_fix, sigma_phi_fix, r_p)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8,
                label=rf'$h_{{\rm LP}} = {hl}~\rm ld$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.axvline(r_p, color='gray', ls=':', lw=1, alpha=0.5)
ax.set_xlabel(r'$r~\rm [ld]$', fontsize=16)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm flat}$', fontsize=16)
ax.set_title(rf'$\rm Vary~h_{{LP}}~(h_p={h_p_fix},~\sigma_r={sigma_r_fix},~\sigma_\phi={sigma_phi_fix})$',
             fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

plt.tight_layout()
plt.savefig('plots/test_pillar_param_scan.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved plots/test_pillar_param_scan.png")

# ======================================================
# Summary: peak enhancement vs shadow depth
# ======================================================
print("\n--- Summary: peak boost at r_p vs shadow depth at r > r_p ---")
print(f"{'h_p':>6s} {'sig_r':>6s} {'sig_phi':>7s} {'h_LP':>6s} {'peak_boost':>11s} {'shadow_min':>11s} {'aspect':>8s}")
for h_p in [0.1, 0.2, 0.5]:
    for sr in [0.2, 0.5]:
        for hl in [0.2, 0.5, 1.0]:
            r, T_ratio, _, _ = get_T_ratio(hl, h_p, sr, 0.2, r_p)
            ir_p = np.argmin(np.abs(r - r_p))
            peak = np.max(T_ratio[max(0,ir_p-5):ir_p+5])
            shadow_region = T_ratio[ir_p+5:]
            shadow_min = np.min(shadow_region) if len(shadow_region) > 0 else 1.0
            aspect = h_p / sr
            print(f"{h_p:6.2f} {sr:6.2f} {0.2:7.2f} {hl:6.2f} {peak:11.4f} {shadow_min:11.4f} {aspect:8.2f}")
