#!/usr/bin/env python3
"""moc_spectrum.py - predicted observed line profiles and spectrum from a
finished MOCASSIN run, compared with the pipeline (Cloudy slabs).

MOCASSIN lines (per-cell luminosities from mocassinPlot, read by
transient.moc_maps.load_run):
  1. Line-of-sight velocity of every sim cell toward the observer at
     i = INCL_DEG (vlos_cube, exactly as the pipeline; positive =
     approaching), resampled onto the run's Cartesian grid as n_H-weighted
     means: n_H v and n_H v^2 go through cmi_convert.resample like n_H,
     so the velocity spread inside a Cartesian cell is kept as a
     broadening sigma_sub. A coarsened run (built from the fine grid with
     the (1/4, 1/2, 1/4) tent average per axis) gets the same coarsening;
     which case applies is detected by matching the run's own grid.npz.
  2. Observer side only (z > 0 plus half the midplane cell), as in
     moc_maps: the far side is hidden behind the opaque midplane.
  3. Each cell's line is a Gaussian at lambda0 (1 - v/c) with
     sigma^2 = sigma_sub^2 + k Te / m_ion, integrated over the bins of
     the pipeline's wavelength grid (~200 km/s), so both spectra share
     one grid. Rest wavelengths follow Cloudy (vacuum below 2000 A, air
     above), as the pipeline's grid does.

Pipeline: compute_rays + per-ray Doppler assembly (as sim_full_spectrum),
split into its reprocessed continuum and line spectra, plus the thermal
(sim T_eff photosphere) and direct-AGN components (thermal_and_agn).

Outputs (in <run>/plots):
  moc_line_profiles.png  L_lambda vs velocity for H alpha, H beta,
                         C IV 1549, Mg II 2798, [O III] 5007: MOCASSIN vs
                         the pipeline's line spectrum (which includes
                         every Cloudy line in the window, e.g. [N II]
                         next to H alpha)
  moc_full_spectrum.png  lambda L_lambda 1000-11000 A: pipeline total,
                         pipeline continuum only, and pipeline continuum
                         + the MOCASSIN lines
and prints line luminosities and FWHMs for both.

Run:  python3 -m transient.moc_spectrum [--run /data2/jhuang/runs/mocassin/moc_coarse]
"""

import argparse
import os

import numpy as np
from scipy.special import erf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import ScalarFormatter, NullFormatter  # noqa: E402

from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, vlos_cube, CONFIG, C_KMS)
from transient.sim_cloudy_spectrum import load_grid, compute_rays  # noqa: E402
from transient.sim_full_spectrum import thermal_and_agn  # noqa: E402
from transient.cmi_convert import resample  # noqa: E402
from transient.moc_convert import grid_half  # noqa: E402
from transient.moc_maps import load_run, DEFAULT_RUN  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

DUMP = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                    'sim', 'disk.out1.00012.athdf')
KB, AMU = 1.380649e-16, 1.66054e-24
# plot.out column -> (rest wavelength [A, Cloudy convention], ion mass [amu])
# for a run with the original seven-line plot.in; runs whose plot.in was
# written by transient.moc_plotin take every line from input/plot_lines.txt
LINE_DATA = {'Ha': (6562.80, 1.008), 'Hb': (4861.33, 1.008),
             'CIV1551': (1550.77, 12.011), 'CIV1548': (1548.19, 12.011),
             'MgII2804': (2802.70, 24.305), 'MgII2796': (2795.53, 24.305),
             'OIII5007': (5006.84, 15.999)}
# profile panels: (title, plot.out columns, velocity zero point [A])
PANELS = [(r'H\alpha', ['Ha'], 6562.80),
          (r'H\beta', ['Hb'], 4861.33),
          (r'C\,IV\,1549', ['CIV1548', 'CIV1551'], 1548.19),
          (r'Mg\,II\,2798', ['MgII2796', 'MgII2804'], 2795.53),
          (r'[O\,III]\,5007', ['OIII5007'], 5006.84)]
VWIN = 10000.0          # km/s half-width of the profile panels
VINT = 6500.0           # km/s half-width for line luminosities (MOCASSIN
                        # |v| < 5600; wider windows pull in other pipeline
                        # lines, e.g. H beta and [O III] 4959 near 5007)


def tent_coarsen(a):
    """(1/4, 1/2, 1/4) tent average by 2 per axis (coarse i on fine 2i),
    edge weights renormalised: the coarsening that built moc_coarse."""
    for ax in range(3):
        a = np.moveaxis(a, ax, 0)
        n = a.shape[0]
        out = np.zeros(((n + 1) // 2,) + a.shape[1:])
        for i in range(out.shape[0]):
            num, wsum = 0.0, 0.0
            for j, w in ((2 * i - 1, .25), (2 * i, .5), (2 * i + 1, .25)):
                if 0 <= j < n:
                    num = num + w * a[j]
                    wsum += w
            out[i] = num / wsum
        a = np.moveaxis(out, 0, ax)
    return a


def velocity_grid(sim, phys, run, shape):
    """n_H-weighted v_los and sub-cell sigma [km/s] on the run's grid."""
    _, (r_o, t_o, p_o) = vlos_cube(sim, CONFIG['INCL_DEG'])
    vlos = phys['vr'] * r_o + phys['vt'] * t_o + phys['vp'] * p_o
    half = grid_half(sim)
    nH_run = np.load(os.path.join(run, 'grid.npz'))['nH']

    def sums(ncell):
        out = [resample(sim, dict(phys, nH=phys['nH'] * vlos ** p), ncell,
                        half, empty=0.0)[1] for p in (0, 1, 2)]
        return out

    for label, ncell, post in (
            ('direct', shape, lambda a: a),
            ('tent-coarsened', [2 * n - 1 for n in shape], tent_coarsen)):
        n0, n1, n2 = (post(a) for a in sums(ncell))
        err = np.max(np.abs(n0 - nH_run)) / nH_run.max()
        if err < 1e-5:
            print(f"velocity grid: {label} resample of {ncell} matches the "
                  f"run's density (max rel. diff {err:.1e})")
            break
    else:
        raise SystemExit("could not reproduce the run's grid.npz density")
    safe = np.where(n0 > 0, n0, 1.0)
    v = np.where(n0 > 0, n1 / safe, 0.0)
    sig = np.sqrt(np.clip(np.where(n0 > 0, n2 / safe, 0.0) - v ** 2, 0, None))
    return v, sig


def binned_lines(wave, L, lam_c, sig_lam):
    """Sum of Gaussians (total L each) integrated over the bins of wave."""
    e = np.concatenate([[1.5 * wave[0] - 0.5 * wave[1]],
                        0.5 * (wave[1:] + wave[:-1]),
                        [1.5 * wave[-1] - 0.5 * wave[-2]]])
    out = np.zeros_like(wave)
    lo, hi = np.searchsorted(e, [lam_c.min() - 6 * sig_lam.max(),
                                 lam_c.max() + 6 * sig_lam.max()])
    lo, hi = max(lo - 1, 0), min(hi + 1, len(e) - 1)
    ee = e[lo:hi + 1]
    cdf = 0.5 * (1 + erf((ee[None, :] - lam_c[:, None])
                         / (np.sqrt(2) * sig_lam[:, None])))
    out[lo:hi] = (L[:, None] * np.diff(cdf, axis=1)).sum(axis=0)
    return out / np.diff(e)                       # erg/s/A


def vac_to_air(lam):
    """Vacuum -> air wavelength [A] (Morton 2000), as Cloudy uses > 2000 A."""
    s2 = (1e4 / lam) ** 2
    n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
    return np.where(lam > 2000, lam / n, lam)


def line_data(d):
    """name -> (rest lambda [A, Cloudy convention], ion mass [amu])."""
    t = d['line_table']
    if 'lam' not in t[0]:
        return LINE_DATA
    return {r['name']: (float(r['lam'] if r['air'] else vac_to_air(r['lam'])),
                        r['mass']) for r in t}


def moc_lines(d, v, sig, wave):
    """Observer-side MOCASSIN line spectra on wave, per plot.out column."""
    z = d['axes'][2]
    zside = np.where(z > 0, 1.0, 0.0)
    zside[np.argmin(np.abs(z))] = 0.5
    sel = d['active'] & (zside[None, None, :] > 0)
    wz = np.broadcast_to(zside[None, None, :], sel.shape)[sel]
    vv, ss, Te = v[sel], sig[sel], d['Te'][sel]
    spec = {}
    for name, (lam0, mass) in line_data(d).items():
        L = d['lines'][name][sel] * 1e36 * wz                 # erg/s
        s_kms = np.sqrt(ss ** 2 + KB * Te / (mass * AMU) / 1e10)
        keep = L > 0
        spec[name] = binned_lines(wave, L[keep],
                                  lam0 * (1 - vv[keep] / C_KMS),
                                  lam0 * s_kms[keep] / C_KMS)
    return spec


def pipeline(sim, phys):
    """Pipeline spectra at INCL_DEG: reprocessed continuum and lines
    (Doppler-assembled per ray as sim_full_spectrum), thermal, AGN."""
    grid = load_grid()
    wave = grid[0]
    rays = compute_rays(sim, phys, CONFIG, grid)
    dA, shift = rays['dA'], 1.0 - rays['v_ray'] / C_KMS
    L_cont, L_line = np.zeros_like(wave), np.zeros_like(wave)
    for i in np.nonzero(dA)[0]:
        ws = wave * shift[i]
        L_cont += dA[i] * np.interp(wave, ws, rays['F_cont'][i] / ws,
                                    left=0, right=0)
        L_line += dA[i] * np.interp(wave, ws, rays['F_line'][i] / ws,
                                    left=0, right=0)
    L_th_base, L_agn = thermal_and_agn(sim, phys, CONFIG, wave)
    L_th = 4.0 * np.cos(np.radians(CONFIG['INCL_DEG'])) * L_th_base
    return wave, dict(cont=L_cont, line=L_line, th=L_th, agn=L_agn)


def peaks(vel, prof):
    """Velocities of the blue and red maxima of a double-peaked profile."""
    b, r = vel < 0, vel > 0
    return vel[b][np.argmax(prof[b])], vel[r][np.argmax(prof[r])]


def ticks(ax):
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default=DEFAULT_RUN)
    ap.add_argument('--dump', default=DUMP)
    a = ap.parse_args()
    d = load_run(a.run)
    if 'lines' not in d:
        raise SystemExit('no output/plot.out: run mocassinPlot first')
    sim = load_sim(a.dump)
    phys = to_physical(sim, CONFIG)
    v, sig = velocity_grid(sim, phys, a.run, d['active'].shape)
    act = d['active']
    print(f"cell v_los: |v| median {np.median(np.abs(v[act])):.0f} km/s, "
          f"max {np.abs(v[act]).max():.0f}; sub-cell sigma median "
          f"{np.median(sig[act]):.0f} km/s")
    wave, P = pipeline(sim, phys)
    M = moc_lines(d, v, sig, wave)
    i_deg = CONFIG['INCL_DEG']

    # ---- profiles ----
    pdir = os.path.join(a.run, 'plots')
    os.makedirs(pdir, exist_ok=True)
    fig, axs = plt.subplots(1, len(PANELS), figsize=(5.2 * len(PANELS), 5))
    print(f"{'line':12s} {'L_MOC':>10s} {'L_pipe':>10s} {'ratio':>6s}   "
          f"peaks MOC / pipe [km/s]   (L in erg/s, |v| < {VINT:.0f})")
    for n, (ax, (title, cols, lam_ref)) in enumerate(zip(axs, PANELS)):
        vel = (wave - lam_ref) / lam_ref * C_KMS            # + = redshift
        win = np.abs(vel) < VWIN
        m = sum(M[c] for c in cols)
        dl = np.gradient(wave)
        wi = np.abs(vel) < VINT
        Lm, Lp = np.sum((m * dl)[wi]), np.sum((P['line'] * dl)[wi])
        pm, pp = peaks(vel[wi], m[wi]), peaks(vel[wi], P['line'][wi])
        name = title.replace('\\,', ' ').replace('\\', '')
        print(f"{name:12s} {Lm:10.3e} {Lp:10.3e} {Lm / Lp:6.2f}   "
              f"{pm[0]:+6.0f} {pm[1]:+6.0f} / {pp[0]:+6.0f} {pp[1]:+6.0f}")
        ax.plot(vel[win], P['line'][win], '-', color='gray', lw=2,
                drawstyle='steps-mid', label=r'$\rm pipeline~(Cloudy)$')
        ax.plot(vel[win], m[win], '-', color='crimson', lw=2,
                drawstyle='steps-mid', label=r'$\rm MOCASSIN$')
        for c in cols[1:]:                      # second doublet component
            dv = (LINE_DATA[c][0] - lam_ref) / lam_ref * C_KMS
            ax.axvline(dv, color='k', ls=':', lw=1)
        ax.axvline(0, color='k', ls=':', lw=1)
        ax.set_title(rf'$\rm {title}$')
        ax.set_xlabel(r'$v~[\rm km~s^{-1}]$')
        if n == 0:
            ax.set_ylabel(r'$L_\lambda~[\rm erg~s^{-1}~\AA^{-1}]$')
            ax.legend(fontsize=13, frameon=False, loc='upper left')
        ax.set_xlim(-VWIN, VWIN)
        ax.set_ylim(0, None)
        ticks(ax)
    fig.suptitle(rf'$\rm line~profiles,~i={i_deg:.0f}^\circ~'
                 rf'(observer~side;~+v=\rm redshift)$')
    fig.tight_layout()
    out = os.path.join(pdir, 'moc_line_profiles.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")

    # ---- full spectrum ----
    L_moc = sum(M.values())
    base = P['cont'] + P['th'] + P['agn']
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    ax.plot(wave, wave * (base + P['line']), '-', color='black', lw=2.5,
            drawstyle='steps-mid', label=r'$\rm pipeline~total$')
    ax.plot(wave, wave * (base + L_moc), '-', color='crimson', lw=2,
            alpha=0.85, drawstyle='steps-mid',
            label=r'$\rm pipeline~continuum + MOCASSIN~lines$')
    # components (colours as sim_full_spectrum)
    ax.plot(wave, wave * (P['cont'] + P['line']), '-', color='darkorange',
            lw=1.8, alpha=0.85, drawstyle='steps-mid',
            label=r'$\rm reprocessed~(Cloudy~continuum + lines)$')
    ax.plot(wave, wave * np.where(L_moc > 0, L_moc, np.nan), '-',
            color='purple', lw=1.8, drawstyle='steps-mid',
            label=r'$\rm MOCASSIN~lines~alone$')
    ax.plot(wave, wave * P['th'], '-', color='seagreen', lw=2.2, alpha=0.9,
            label=r'${\rm thermal~(sim~}T_{\rm eff}{\rm ,~disk+TDE+arms)}$')
    ax.plot(wave, wave * P['agn'], '--', color='royalblue', lw=2.2,
            alpha=0.9, label=r'$\rm AGN~direct~(Cloudy~incident~SED)$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1000, 11000)
    msel = (wave >= 1000) & (wave <= 11000)
    ymax = (wave * (base + np.maximum(P['line'], L_moc)))[msel].max()
    ax.set_ylim(ymax / 3e3, ymax * 3)
    ax.set_xticks([1000, 2000, 3000, 5000, 10000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'$\rm wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\lambda L_\lambda~[\rm erg~s^{-1}]$', fontsize=17)
    ax.set_title(rf'$i={i_deg:.0f}^\circ;~\rm MOCASSIN:~{len(M)}~lines~'
                 r'(no~Ly\alpha,~no~Fe)$', fontsize=16)
    ax.legend(fontsize=12, frameon=False, loc='upper right')
    ticks(ax)
    out = os.path.join(pdir, 'moc_full_spectrum.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == '__main__':
    main()
