#!/usr/bin/env python3
"""
fit_T2_profile.py - fully NONLINEAR arbitrary flare temperature profile:
the flare-state disk has a free positive temperature profile T2(r) (no
perturbation bound), and the exact Planck-difference SED

    dF(lam) = sum_g W_g [ B(lam, T2_g) - B(lam, Tq_g) ]

(radial groups g, weights and quiescent temperatures Tq from the model
disk) is fitted to the observed 2019-2000 difference, keeping the
adopted raised-arm occultation. Parameters: ln T2 per group (40),
optimized with scipy least_squares (analytic-free-form, mild ln-space
smoothness, multiple starts).

Outputs: printed chi2 per smoothness/start; three-panel spectrum figure
(transient_T2_model.png) and temperature-profile figure
(transient_T2_profile.png).

Run:  python3 transient/fit_T2_profile.py
"""

import copy
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

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
NGROUP = 40
SMOOTH = [1.0, 5.0, 20.0]      # ln-space curvature penalty weights
LNT_LO, LNT_HI = np.log(300.0), np.log(5e5)
NR, NPHI = 200, 400


def main():
    # ---- data, geometry, per-ring weights ----
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

    edges = np.logspace(np.log10(rr[0]), np.log10(rr[-1]), NGROUP + 1)
    gidx = np.clip(np.digitize(rr, edges) - 1, 0, NGROUP - 1)
    used = np.array([np.any(gidx == j) for j in range(NGROUP)])
    W_g = np.array([W_r[gidx == j].sum() for j in range(NGROUP)])[used]
    Tq_g = np.array([np.average(T_r[gidx == j], weights=W_r[gidx == j])
                     if np.any(gidx == j) else np.nan
                     for j in range(NGROUP)])[used]
    r_g = np.sqrt(edges[:-1] * edges[1:])[used]
    ng = used.sum()

    # Planck matrix helpers (B evaluated per wavelength for a T vector)
    def Bmat(Tvec):
        return np.array([quiet.planck_function(w, Tvec) for w in WL])

    B_q = Bmat(Tq_g)                                  # (nwl, ng)
    F_base = (B_q * W_g[None, :]).sum(axis=1)

    wgt = 1.0 / derr ** 2
    sqw = np.sqrt(wgt[good])
    ndof = good.sum()
    eps = 1e-3

    def model_diff(lnT):
        T2 = np.exp(lnT)
        F2 = (Bmat(T2) * W_g[None, :]).sum(axis=1)
        return (F2 - F_base) * fnorm - occ_def

    def make_fun(lam_s):
        def fun(lnT):
            res = (model_diff(lnT)[good] - ddiff[good]) * sqw
            smooth = lam_s * np.diff(lnT, 2)
            return np.concatenate([res, smooth])

        def jac(lnT):
            T2 = np.exp(lnT)
            dB = (Bmat(T2 * (1 + eps)) - Bmat(T2)) / eps   # dB/dlnT
            Jd = (dB * W_g[None, :])[good] * fnorm * sqw[:, None]
            Js = lam_s * np.diff(np.eye(ng), 2, axis=0)
            return np.vstack([Jd, Js])
        return fun, jac

    starts = {
        'quiescent': np.log(Tq_g),
        'gauss': np.log(Tq_g * (1 + 0.5 * np.exp(-0.5 * ((r_g - 2.) / 1.)
                                                 ** 2))),
    }
    print(f"  {'smooth':>7} {'start':>10} {'chi2/dof':>9}")
    results = []
    for lam_s in SMOOTH:
        fun, jac = make_fun(lam_s)
        for sname, p0 in starts.items():
            sol = least_squares(fun, p0, jac=jac, method='trf',
                                bounds=(LNT_LO, LNT_HI), max_nfev=400)
            chi2 = float(np.sum(
                ((model_diff(sol.x)[good] - ddiff[good]) ** 2 * wgt[good])))
            print(f"  {lam_s:>7.1f} {sname:>10} {chi2/ndof:>9.1f}",
                  flush=True)
            results.append((chi2, lam_s, sname, sol.x))

    chi2, lam_s, sname, lnT_best = min(results, key=lambda r: r[0])
    T2_best = np.exp(lnT_best)
    dF_best = model_diff(lnT_best) + occ_def
    model_best = dF_best - occ_def
    print(f"\n  best: smooth = {lam_s:g}, start = {sname}, "
          f"chi2/dof = {chi2/ndof:.1f} "
          f"(adopted Balmer-edge model: 31.6)")

    # ---- UV-upweighted free fit (x10 chi2 weight blueward of 3050 A) ----
    wgt_uv = wgt * np.where(WL < 3050., 10.0, 1.0)
    sqw_uv = np.sqrt(wgt_uv[good])

    def make_fun_uv(lam_s):
        def fun(lnT):
            res = (model_diff(lnT)[good] - ddiff[good]) * sqw_uv
            return np.concatenate([res, lam_s * np.diff(lnT, 2)])

        def jac(lnT):
            T2 = np.exp(lnT)
            dB = (Bmat(T2 * (1 + eps)) - Bmat(T2)) / eps
            Jd = (dB * W_g[None, :])[good] * fnorm * sqw_uv[:, None]
            Js = lam_s * np.diff(np.eye(ng), 2, axis=0)
            return np.vstack([Jd, Js])
        return fun, jac

    fun_uv, jac_uv = make_fun_uv(1.0)
    best_uv = None
    for sname, p0 in starts.items():
        sol = least_squares(fun_uv, p0, jac=jac_uv, method='trf',
                            bounds=(LNT_LO, LNT_HI), max_nfev=400)
        obj = float(np.sum(sol.fun ** 2))
        if best_uv is None or obj < best_uv[0]:
            best_uv = (obj, sol.x)
    lnT_uv = best_uv[1]
    T2_uv = np.exp(lnT_uv)
    model_uv = model_diff(lnT_uv)
    dF_uv = model_uv + occ_def
    chi2_uvfit = float(np.sum((model_uv[good] - ddiff[good]) ** 2
                              * wgt[good]))
    print(f"  UV-weighted free T2: plain chi2/dof = {chi2_uvfit/ndof:.1f}")

    # Gaussian-bump temperature model as a fixed reference
    T2_gref = T_r * (1.0 + 0.5 * np.exp(-0.5 * ((rr - 2.0) / 1.0) ** 2))
    dF_gref = np.array([np.sum(W_r * (quiet.planck_function(w, T2_gref)
                                      - quiet.planck_function(w, T_r)))
                        for w in WL]) * fnorm
    chi2_g = float(np.sum(((dF_gref - occ_def)[good] - ddiff[good]) ** 2
                          * wgt[good]))
    print(f"  Gaussian-bump reference (r0=2, a=0.5): "
          f"chi2/dof = {chi2_g/ndof:.1f}")

    # SINGLE Gaussian island: T2(r) IS one Gaussian (not on top of the
    # old profile; the rest of the disk is cold). (T_pk, r0, sigma) fitted.
    def island(p):
        tpk, r0i, sgi = np.exp(p[0]), p[1], np.exp(p[2])
        return np.maximum(tpk * np.exp(-0.5 * ((rr - r0i) / sgi) ** 2),
                          300.0)

    def dflux_fine(T2):
        return np.array([np.sum(W_r * (quiet.planck_function(w, T2)
                                       - quiet.planck_function(w, T_r)))
                         for w in WL]) * fnorm

    def resid_isl(p):
        return (dflux_fine(island(p)) - occ_def - ddiff)[good] * sqw_g

    sqw_g = np.sqrt(wgt[good])
    best_isl = None
    for r0i in (1.0, 2.0, 4.0):
        sol = least_squares(resid_isl, [np.log(18000.), r0i, np.log(1.5)],
                            bounds=([np.log(2000.), 0.05, np.log(0.1)],
                                    [np.log(2e5), 30.0, np.log(20.0)]),
                            diff_step=1e-3, max_nfev=200)
        chi2_i = float(np.sum(sol.fun ** 2))
        if best_isl is None or chi2_i < best_isl[0]:
            best_isl = (chi2_i, sol.x)
    chi2_i, p_isl = best_isl
    T2_isl = island(p_isl)
    dF_isl = dflux_fine(T2_isl)
    print(f"  single Gaussian island: T_pk = {np.exp(p_isl[0])/1e3:.1f} kK,"
          f" r0 = {p_isl[1]:.2f} ld, sigma = {np.exp(p_isl[2]):.2f} ld,"
          f" chi2/dof = {chi2_i/ndof:.1f}")
    isl_lab_top = (rf'$\rm flare~(single~Gaussian~island,~'
                   rf'{np.exp(p_isl[0])/1e3:.0f}~kK~at~'
                   rf'{p_isl[1]:.0f}~ld)$')
    isl_lab = r'$\rm single~Gaussian~island$'

    # UV-upweighted island fit (same x10 blue weighting as the free fit)
    def resid_isl_uv(p):
        return (dflux_fine(island(p)) - occ_def - ddiff)[good] * sqw_uv

    best_isl_uv = None
    for r0i in (1.0, 2.0, 4.0):
        sol = least_squares(resid_isl_uv,
                            [np.log(18000.), r0i, np.log(1.5)],
                            bounds=([np.log(2000.), 0.05, np.log(0.1)],
                                    [np.log(2e5), 30.0, np.log(20.0)]),
                            diff_step=1e-3, max_nfev=200)
        obj = float(np.sum(sol.fun ** 2))
        if best_isl_uv is None or obj < best_isl_uv[0]:
            best_isl_uv = (obj, sol.x)
    p_isl_uv = best_isl_uv[1]
    T2_isl_uv = island(p_isl_uv)
    dF_isl_uv = dflux_fine(T2_isl_uv)
    chi2_iuv = float(np.sum(((dF_isl_uv - occ_def)[good] - ddiff[good])
                            ** 2 * wgt[good]))
    print(f"  island, UV-weighted: T_pk = "
          f"{np.exp(p_isl_uv[0])/1e3:.1f} kK, r0 = {p_isl_uv[1]:.2f} ld,"
          f" sigma = {np.exp(p_isl_uv[2]):.2f} ld, plain chi2/dof = "
          f"{chi2_iuv/ndof:.1f}")

    # ---- adopted reference ----
    _, _, _, _, occ_def_s, A, flare_ref = tr.setup()
    shape = np.array([flare_ref._arm_emission_shape(
        w, ANA['te'], ANA['fbb'], ANA['vs'], ANA['pas']) for w in WL])
    E = A * shape
    num = np.sum((ddiff[good] + occ_def_s[good]) * E[good] * wgt[good])
    den = np.sum(E[good] ** 2 * wgt[good])
    model_edge = max(0.0, num / den) * E - occ_def_s

    # ---- UV-side metrics (the discriminating region) ----
    blue = good & (WL < 3050.)
    nblue = blue.sum()
    i2550 = np.argmin(np.abs(WL - 2550.))
    print(f"\n  UV side (< 3050 A, n = {nblue}); data dF(2550) = "
          f"{ddiff[i2550]:+.1f}")
    for lab, m in (('Balmer continuum', model_edge),
                   ('free T2(r)', model_best),
                   ('free T2, UV-weighted', model_uv),
                   ('Gaussian bump', dF_gref - occ_def),
                   ('single island', dF_isl - occ_def),
                   ('island, UV-weighted', dF_isl_uv - occ_def)):
        c2b = float(np.sum((m[blue] - ddiff[blue]) ** 2 * wgt[blue]))
        print(f"  {lab:>17}: chi2_UV/n = {c2b/nblue:>8.1f}, "
              f"model dF(2550) = {m[i2550]:+.1f}")

    flam_quiet = q * scale
    flam_T2 = occ0 + dF_best
    flam_edge = occ0 + (model_edge + occ_def_s)

    # ---- figure 1: three-panel spectra ----
    binned = {}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        ep = load_epoch(name)
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        binned[name] = (wb / (1.0 + tr.Z), fb, eb)
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
    ax1.plot(WL, flam_quiet, '-', color='darkorange', lw=3, alpha=0.9,
             label=r'$\rm model~quiescent~(no~arm)$')
    ax1.plot(WL, flam_edge, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm flare~(include~Balmer~continuum)$')
    ax1.plot(WL, occ0 + dF_gref, '-', color='purple', lw=3, alpha=0.85,
             label=r'$\rm flare~(Gaussian~bump,~r_0=2~ld,~a=0.5)$')
    ax1.plot(WL, occ0 + dF_uv, '-', color='deeppink', lw=3, alpha=0.9,
             label=r'$\rm flare~(free~T_2(r))$')
    ax1.plot(WL, occ0 + dF_isl_uv, '-', color='darkturquoise', lw=3,
             alpha=0.9,
             label=r'$\rm flare~(single~Gaussian)$')

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
             label=r'$\rm Balmer~continuum$')
    ax2.plot(WL, dF_gref - occ_def, '-', color='purple', lw=3,
             alpha=0.85, label=r'$\rm Gaussian~bump$')
    ax2.plot(WL, model_uv, '-', color='deeppink', lw=3, alpha=0.9,
             label=r'$\rm free~T_2(r)$')
    ax2.plot(WL, dF_isl_uv - occ_def, '-', color='darkturquoise', lw=3,
             alpha=0.9, label=r'$\rm single~Gaussian$')

    r_val = f2019 / f2000b
    r_err = r_val * np.sqrt((e2019 / f2019) ** 2 + (e2000 / f2000b) ** 2)
    ax3.fill_between(grid_rest, r_val - r_err, r_val + r_err, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax3.plot(grid_rest, r_val, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019/2000~(data)$')
    ax3.plot(WL, flam_edge / flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm Balmer~continuum$')
    ax3.plot(WL, (occ0 + dF_gref) / flam_quiet, '-', color='purple',
             lw=2, alpha=0.85, label=r'$\rm Gaussian~bump$')
    ax3.plot(WL, (occ0 + dF_uv) / flam_quiet, '-', color='deeppink',
             lw=3, alpha=0.9, label=r'$\rm free~T_2(r)$')
    ax3.plot(WL, (occ0 + dF_isl_uv) / flam_quiet, '-',
             color='darkturquoise', lw=3, alpha=0.9,
             label=r'$\rm single~Gaussian$')

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
    out = os.path.join(PLOTDIR, 'transient_T2_model.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- figure 2: the temperature profiles ----
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(rr, T_r, '-', color='darkorange', lw=3, alpha=0.9,
            label=r'$\rm quiescent~T(r)$')
    ax.plot(rr, T2_gref, '-', color='purple', lw=3, alpha=0.85,
            label=r'$\rm Gaussian~bump~(r_0=2~ld,~a=0.5)$')
    ax.plot(r_g, T2_uv, '-', color='deeppink', lw=3, alpha=0.9,
            label=r'$\rm free~T_2(r)$')
    ax.plot(rr, T2_isl_uv, '-', color='darkturquoise', lw=3, alpha=0.9,
            label=r'$\rm single~Gaussian$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(rr[0], 30)
    ax.set_ylim(1e3, 3e5)
    ax.legend(fontsize=14, frameon=False)
    ax.set_xlabel(r'$\rm radius~[light~days]$', fontsize=18)
    ax.set_ylabel(r'$\rm T~[K]$', fontsize=18)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out2 = os.path.join(PLOTDIR, 'transient_T2_profile.png')
    plt.savefig(out2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out2}")


if __name__ == '__main__':
    main()
