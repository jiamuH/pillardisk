#!/usr/bin/env python3
"""
fit_dT_profile.py - OPTIMAL temperature-change profile: allow an
arbitrary (positive or negative) perturbation dT(r) of the quiet disk's
temperature profile and solve linearly for the dT(r) that best fits the
observed 2019-2000 difference (with the adopted raised-arm occultation
kept). In the perturbative limit
    dF(lam) = sum_r W(r) dB/dT(lam, T(r)) dT(r),
so dT(r) follows from ridge-regularized weighted least squares; the
result is then re-evaluated with the exact Planck difference.

Outputs:
  - printed chi2 vs regularization strength (linear + exact)
  - three-panel spectrum figure (data + model quiescent + dT-profile
    flare model + adopted Balmer-edge model): transient_dT_model.png
  - the profile itself (dT/T vs r): transient_dT_profile.png

Run:  python3 transient/fit_dT_profile.py
"""

import copy
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import transient.tune_redside as tr  # noqa: E402
from transient.transient_disk import disks_from_config  # noqa: E402
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, LINE_WINDOWS, PLOTDIR)
from pillar_disk import C, DAY, FNU_TO_MJY  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

ANA = dict(te=12000., fbb=0.10, vs=55000., pas=0.0)   # adopted reference
NGROUP = 40            # radial groups (log-spaced) for dT(r)
BOUNDS = [0.2, 0.5, 1.0]   # allowed |dT/T| per radial group
ALPHA_S = 3.0              # smoothness weight (relative)
NR, NPHI = 200, 400


def main():
    # ---- data, geometry, per-ring weights (as in fit_ring_heating) ----
    WL, ddiff, derr, f2000, good = tr.data_on_grid()
    nw = tr.CFG['transient']['normalization']['rest_window']
    normsel = (WL >= nw[0]) & (WL <= nw[1])
    c0 = copy.deepcopy(tr.CFG)
    c0['transient']['pillar']['pillar_temp'] = 0.0
    flare0, quiet, _ = disks_from_config(c0, nr=NR, nphi=NPHI)
    q = quiet.compute_sed(WL) * tr.TO_FLAM
    scale = np.median(f2000[normsel]) / np.median(q[normsel])
    occ0 = flare0.compute_sed(WL) * tr.TO_FLAM * scale
    occ_def = q * scale - occ0

    weight, _ = quiet._sed_weights(False)
    W_r = weight.sum(axis=1)
    r_2d, phi_2d = np.meshgrid(quiet.r, quiet.phi, indexing='ij')
    T_r = quiet.get_temperature(r_2d, phi_2d)[:-1, :].mean(axis=1)
    rr = quiet.r[:-1]
    ld_to_cm = C * DAY
    fnorm = (ld_to_cm ** 2 / (quiet.d * ld_to_cm) ** 2 * FNU_TO_MJY
             * tr.TO_FLAM * scale)

    # ---- radial groups and dB/dT basis ----
    edges = np.logspace(np.log10(rr[0]), np.log10(rr[-1]), NGROUP + 1)
    gidx = np.clip(np.digitize(rr, edges) - 1, 0, NGROUP - 1)
    K = np.zeros((len(WL), NGROUP))
    eps = 1e-3
    for iw, w in enumerate(WL):
        dBdT = (quiet.planck_function(w, T_r * (1 + eps))
                - quiet.planck_function(w, T_r)) / (eps * T_r)
        np.add.at(K[iw], gidx, W_r * dBdT)
    K *= fnorm
    T_g = np.array([np.median(T_r[gidx == j]) if np.any(gidx == j)
                    else np.nan for j in range(NGROUP)])
    r_g = np.sqrt(edges[:-1] * edges[1:])
    used = np.array([np.any(gidx == j) for j in range(NGROUP)])
    K = K[:, used]
    T_g, r_g = T_g[used], r_g[used]
    ng = used.sum()

    # ---- bounded solve: model = K' u - occ_def, u = dT/T per group ----
    from scipy.optimize import lsq_linear
    wgt = 1.0 / derr ** 2
    target = ddiff + occ_def
    Kp = K * T_g[None, :]                           # columns now per dT/T
    Kw = Kp[good] * np.sqrt(wgt[good])[:, None]
    tw = target[good] * np.sqrt(wgt[good])
    L = np.diff(np.eye(ng), 2, axis=0)              # 2nd-difference
    srow = np.sqrt(ALPHA_S * np.trace(Kw.T @ Kw) / ng)
    A_des = np.vstack([Kw, srow * L])
    b_des = np.concatenate([tw, np.zeros(L.shape[0])])
    ndof = good.sum()
    print(f"  {'|dT/T|<':>8} {'chi2/dof lin':>13} {'chi2/dof exact':>15} "
          f"{'max|dT/T|':>10}")
    results = []
    for bnd in BOUNDS:
        sol = lsq_linear(A_des, b_des, bounds=(-bnd, bnd))
        u = sol.x
        model_lin = Kp @ u - occ_def
        chi2_lin = float(np.sum((model_lin[good] - ddiff[good]) ** 2
                                * wgt[good]))
        # exact nonlinear evaluation of the same dT(r)/T
        ufull = np.zeros(NGROUP)
        ufull[used] = u
        T2 = np.maximum(T_r * (1.0 + ufull[gidx]), 1000.0)
        dF = np.array([np.sum(W_r * (quiet.planck_function(w, T2)
                                     - quiet.planck_function(w, T_r)))
                       for w in WL]) * fnorm
        model_ex = dF - occ_def
        chi2_ex = float(np.sum((model_ex[good] - ddiff[good]) ** 2
                               * wgt[good]))
        print(f"  {bnd:>8.1f} {chi2_lin/ndof:>13.1f} "
              f"{chi2_ex/ndof:>15.1f} {np.max(np.abs(u)):>10.2f}")
        results.append((chi2_ex, bnd, u * T_g, model_ex, dF))

    best = min(results, key=lambda r: r[0])
    chi2_ex, bnd, c, model_ex, dF = best
    print(f"\n  adopted profile: |dT/T| < {bnd:.1f}, "
          f"chi2/dof (exact) = {chi2_ex/ndof:.1f}")

    # ---- comparison temperature models (exact SEDs) ----
    T_at2 = float(np.interp(2.0, rr, T_r))
    comparisons = {}
    T2_ss = 1.5 * T_r
    comparisons['scaled'] = (T2_ss, 'teal',
                             r'$\rm 1.5\,T(r)~(scaled~SS~disk)$')
    T2_g = T_r * (1.0 + 0.5 * np.exp(-0.5 * ((rr - 2.0) / 1.0) ** 2))
    comparisons['gauss'] = (T2_g, 'chocolate',
                            r'$\rm Gaussian~bump~(r_0=2~ld,~a=0.5)$')
    T2_p = np.where((rr >= 0.7) & (rr <= 3.0), 1.5 * T_at2, T_r)
    comparisons['plateau'] = (T2_p, 'purple',
                              r'$\rm plateau~T=1.5\,T(2~ld),~0.7\mbox{-}3~ld$')
    comp_dF = {}
    for key, (T2c, col, lab) in comparisons.items():
        dFc = np.array([np.sum(W_r * (quiet.planck_function(w, T2c)
                                      - quiet.planck_function(w, T_r)))
                        for w in WL]) * fnorm
        chi2c = float(np.sum(((dFc - occ_def)[good] - ddiff[good]) ** 2
                             * wgt[good]))
        comp_dF[key] = dFc
        print(f"  {key:>8}: chi2/dof (exact) = {chi2c/ndof:.1f}")

    # ---- adopted Balmer-edge reference on WL ----
    _, _, _, _, occ_def_s, A, flare_ref = tr.setup()
    shape = np.array([flare_ref._arm_emission_shape(
        w, ANA['te'], ANA['fbb'], ANA['vs'], ANA['pas']) for w in WL])
    E = A * shape
    num = np.sum((ddiff[good] + occ_def_s[good]) * E[good] * wgt[good])
    den = np.sum(E[good] ** 2 * wgt[good])
    model_edge = max(0.0, num / den) * E - occ_def_s

    # ---- figure 1: three-panel spectra ----
    binned = {}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        ep = load_epoch(name)
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        binned[name] = (wb / (1.0 + tr.Z), fb, eb)
    flam_quiet = q * scale
    flam_dT = occ0 + dF                       # occulted disk + heating
    flam_edge_full = occ0 + (model_edge + occ_def_s)

    fig, axes = plt.subplots(3, 1, figsize=(12, 13), sharex=True,
                             gridspec_kw={'height_ratios': [2.2, 1, 1],
                                          'hspace': 0.06})
    ax1, ax2, ax3 = axes
    colors = {'sdss2000': 'crimson', 'eboss2019': 'black',
              'desi2021': 'royalblue'}
    labels = {'sdss2000': r'$\rm SDSS\mbox{-}2000$',
              'eboss2019': r'$\rm SDSS\mbox{-}2019$',
              'desi2021': r'$\rm DESI\mbox{-}2021$'}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        if name in binned:
            wrb, frb, erb = binned[name]
            ax1.fill_between(wrb, frb - erb, frb + erb, step='mid',
                             color=colors[name], alpha=0.22, lw=0)
            ax1.plot(wrb, frb, drawstyle='steps-mid', color=colors[name],
                     lw=1.3, alpha=0.9, label=labels[name])
    ax1.plot(WL, flam_quiet, '--', color='darkorange', lw=3, alpha=0.9,
             label=r'$\rm model~quiescent~(no~arm)$')
    ax1.plot(WL, flam_edge_full, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm adopted~flare~(Balmer\mbox{-}edge~arm)$')
    ax1.plot(WL, flam_dT, '-', color='mediumvioletred', lw=3, alpha=0.9,
             label=r'$\rm flare~(optimal~\delta T(r))$')
    for key, (T2c, col, lab) in comparisons.items():
        ax1.plot(WL, occ0 + comp_dF[key], '-', color=col, lw=2, alpha=0.85,
                 label=lab)

    grid_rest = binned['sdss2000'][0]
    f2000b = binned['sdss2000'][1]
    e2000 = binned['sdss2000'][2]
    f2019 = np.interp(grid_rest, binned['eboss2019'][0],
                      binned['eboss2019'][1])
    e2019 = np.interp(grid_rest, binned['eboss2019'][0],
                      binned['eboss2019'][2])
    d_err = np.sqrt(e2019 ** 2 + e2000 ** 2)
    ax2.fill_between(grid_rest, f2019 - f2000b - d_err,
                     f2019 - f2000b + d_err, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax2.plot(grid_rest, f2019 - f2000b, drawstyle='steps-mid',
             color='darkgreen', lw=1.5, alpha=0.9,
             label=r'$\rm 2019-2000~(data)$')
    ax2.plot(WL, model_edge, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm adopted$')
    ax2.plot(WL, model_ex, '-', color='mediumvioletred', lw=3, alpha=0.9,
             label=r'$\rm optimal~\delta T(r)$')
    for key, (T2c, col, lab) in comparisons.items():
        ax2.plot(WL, comp_dF[key] - occ_def, '-', color=col, lw=2,
                 alpha=0.85)

    r_val = f2019 / f2000b
    r_err = r_val * np.sqrt((e2019 / f2019) ** 2 + (e2000 / f2000b) ** 2)
    ax3.fill_between(grid_rest, r_val - r_err, r_val + r_err, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax3.plot(grid_rest, r_val, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019/2000~(data)$')
    ax3.plot(WL, flam_edge_full / flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm adopted$')
    ax3.plot(WL, flam_dT / flam_quiet, '-', color='mediumvioletred', lw=3,
             alpha=0.9, label=r'$\rm optimal~\delta T(r)$')
    for key, (T2c, col, lab) in comparisons.items():
        ax3.plot(WL, (occ0 + comp_dF[key]) / flam_quiet, '-', color=col,
                 lw=2, alpha=0.85)

    ax1.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                   fontsize=18)
    ax1.legend(fontsize=12, frameon=False, ncol=2, loc='upper right')
    ax2.axhline(0, color='gray', lw=1)
    ax2.set_ylabel(r'$\Delta f_\lambda$', fontsize=18)
    ax2.legend(fontsize=12, frameon=False, loc='upper right')
    ax3.axhline(1, color='gray', lw=1)
    ax3.set_ylabel(r'$\rm ratio$', fontsize=18)
    ax3.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax3.legend(fontsize=12, frameon=False, loc='upper left')
    for ax in axes:
        for w1_, w2_ in LINE_WINDOWS:
            ax.axvspan(w1_, w2_, color='gray', alpha=0.12, lw=0)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
    ax1.set_xlim(2450, 6000)
    ax1.set_ylim(0, None)
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'transient_dT_model.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- figure 2: the temperature profiles ----
    cfull = np.zeros(NGROUP)
    cfull[used] = c
    T2_opt = np.maximum(T_r + cfull[gidx], 1000.0)
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(rr, T_r, '--', color='darkorange', lw=3, alpha=0.9,
            label=r'$\rm quiescent~T(r)$')
    ax.plot(rr, T2_opt, '-', color='mediumvioletred', lw=3, alpha=0.9,
            label=r'$\rm optimal~\delta T(r)$')
    for key, (T2c, col, lab) in comparisons.items():
        ax.plot(rr, T2c, '-', color=col, lw=2, alpha=0.85, label=lab)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(rr[0], 30)
    ax.set_ylim(1e3, 3e5)
    ax.legend(fontsize=13, frameon=False)
    ax.set_xlabel(r'$\rm radius~[light~days]$', fontsize=18)
    ax.set_ylabel(r'$\rm T~[K]$', fontsize=18)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out2 = os.path.join(PLOTDIR, 'transient_dT_profile.png')
    plt.savefig(out2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out2}")


if __name__ == '__main__':
    main()
