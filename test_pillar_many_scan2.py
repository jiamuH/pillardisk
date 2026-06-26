"""
Fine scan around the best config (h_p=0.5, sigma_r=0.5, sigma_phi=0.3).
Explore nearby parameter space to find optimal enhancement.
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


def add_pillars(disk, N, h_pillar, sigma_r, sigma_phi, r_mean=10.0, sig_r_dist=5.0, seed=42):
    rng = np.random.default_rng(seed)
    r_vals = rng.normal(r_mean, sig_r_dist, N)
    r_vals = np.clip(r_vals, disk.rin + 0.5, disk.rout - 0.5)
    phi_vals = rng.uniform(0, 2*np.pi, N)
    for i in range(N):
        disk.add_pillar(r_pillar=r_vals[i], phi_pillar=phi_vals[i],
                        height=h_pillar, sigma_r=sigma_r, sigma_phi=sigma_phi,
                        modify_height=True, modify_temp=True)


def get_T_ratio(hlamp, h_p, sr, sp):
    disk = build_disk(hlamp=hlamp)
    add_pillars(disk, 100, h_p, sr, sp)
    r, T_with, T_without, _, _, _ = disk.get_temperature_profile()
    return r, T_with / np.maximum(T_without, 1.0)


def style_ax(ax):
    ax.tick_params(direction='in', which='major', length=8, width=1.8, labelsize=13,
                   top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


hlamp = 0.5

# Baseline best: h_p=0.5, sr=0.5, sp=0.3
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Top left: vary h_p around 0.5 (fix sr=0.5, sp=0.3)
ax = axes[0, 0]
sr_fix, sp_fix = 0.5, 0.3
for h_p in [0.3, 0.4, 0.5, 0.7, 1.0]:
    r, T_ratio = get_T_ratio(hlamp, h_p, sr_fix, sp_fix)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8, label=rf'$h_p = {h_p}$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm bowl}$', fontsize=16)
ax.set_title(rf'$\rm Vary~h_p~(\sigma_r={sr_fix},~\sigma_\phi={sp_fix})$', fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

# Top right: vary sr around 0.5 (fix h_p=0.5, sp=0.3)
ax = axes[0, 1]
hp_fix, sp_fix = 0.5, 0.3
for sr in [0.3, 0.4, 0.5, 0.7, 1.0]:
    r, T_ratio = get_T_ratio(hlamp, hp_fix, sr, sp_fix)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8, label=rf'$\sigma_r = {sr}$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm bowl}$', fontsize=16)
ax.set_title(rf'$\rm Vary~\sigma_r~(h_p={hp_fix},~\sigma_\phi={sp_fix})$', fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

# Bottom left: vary sp around 0.3 (fix h_p=0.5, sr=0.5)
ax = axes[1, 0]
hp_fix, sr_fix = 0.5, 0.5
for sp in [0.1, 0.2, 0.3, 0.5, 0.8]:
    r, T_ratio = get_T_ratio(hlamp, hp_fix, sr_fix, sp)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8, label=rf'$\sigma_\phi = {sp}$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.set_xlabel(r'$r~\rm [ld]$', fontsize=16)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm bowl}$', fontsize=16)
ax.set_title(rf'$\rm Vary~\sigma_\phi~(h_p={hp_fix},~\sigma_r={sr_fix})$', fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

# Bottom right: vary h_lamp (fix h_p=0.5, sr=0.5, sp=0.3)
ax = axes[1, 1]
hp_fix, sr_fix, sp_fix = 0.5, 0.5, 0.3
for hl in [0.2, 0.3, 0.5, 0.7, 1.0]:
    r, T_ratio = get_T_ratio(hl, hp_fix, sr_fix, sp_fix)
    ax.semilogx(r, T_ratio, lw=2.5, alpha=0.8, label=rf'$h_{{\rm LP}} = {hl}$')
ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.set_xlabel(r'$r~\rm [ld]$', fontsize=16)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm bowl}$', fontsize=16)
ax.set_title(rf'$\rm Vary~h_{{LP}}~(h_p={hp_fix},~\sigma_r={sr_fix},~\sigma_\phi={sp_fix})$', fontsize=14, pad=10)
ax.legend(fontsize=11)
ax.set_xlim(1, 20)
style_ax(ax)

plt.tight_layout()
plt.savefig('plots/test_pillar_many_scan2.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved plots/test_pillar_many_scan2.png")

# --- Summary table ---
print(f"\n--- Fine scan: mean T_ratio over r=5-15 ld (100 pillars, h_LP={hlamp}) ---")
print(f"{'h_p':>6s} {'sig_r':>6s} {'sig_phi':>7s} {'h_LP':>6s} {'mean_boost':>11s} {'max_boost':>10s} {'min_ratio':>10s}")
for h_p in [0.3, 0.5, 0.7, 1.0]:
    for sr in [0.3, 0.5, 0.7]:
        for sp in [0.2, 0.3, 0.5]:
            r, T_ratio = get_T_ratio(hlamp, h_p, sr, sp)
            mask = (r > 5) & (r < 15)
            mean_b = np.mean(T_ratio[mask])
            print(f"{h_p:6.1f} {sr:6.1f} {sp:7.1f} {hlamp:6.1f} {mean_b:11.4f} {np.max(T_ratio):10.4f} {np.min(T_ratio[r>2]):10.4f}")
