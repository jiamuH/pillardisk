#!/usr/bin/env python3
"""opacity_tables.py - Rosseland mean opacity against temperature at fixed
density from the published tables, for the disk-atmosphere opacity
question (dust below T_sub ~ 1500 K; the hydrogen-ionization rise near
1e4 K, as in the dwarf-nova and Cepheid instabilities).

Tables (X = 0.70, Z = 0.02, Grevesse & Noels 1993 solar mix), downloaded to
TABLE_DIR:
  log T <= T_JOIN: Ferguson et al. 2005, ApJ 623, 585 (molecules, grains),
                   f05/x/g7.02.tron from the Wichita State opacity page
  log T >  T_JOIN: OPAL (Iglesias & Rogers 1996, ApJ 464, 943), GN93hz
                   table 73, from A. Boothroyd's CITA mirror
Both give log kappa_R [cm^2/g] on (log T, log R), R = rho / T6^3. At a
fixed rho, log R = log rho - 3 (log T - 6); below the tables' log R = -8
(hot, dilute: electron scattering) the opacity is held at its log R = -8
value (drawn dotted).

Run:  python3 -m transient.opacity_tables   (writes PLOT)
      python3 -m transient.opacity_tables components   (PLOT_COMP: the
      table total against analytic LTE gas components - H- bound-free and
      free-free, H bound-free and free-free, electron scattering - as
      additive shares of the analytic kappa_R; Saha with H, He and lumped
      metal donors; no molecules, dust or lines)
"""

import os

import numpy as np

TABLE_DIR = '/data2/jhuang/data/opacity_tables'
F05 = os.path.join(TABLE_DIR, 'f05', 'x', 'g7.02.tron')
OPAL = os.path.join(TABLE_DIR, 'opal', 'GN93hz')
OPAL_TABLE = 73                         # X = 0.70, Z = 0.02
# the main plot's metallicities: (Z, F05 file, OPAL table, line style,
# label); X = 0.70 for both (5 Z_sun: Y = 0.20)
METALS = ((0.02, 'g7.02.tron', 73, '-', r'$\bm{Z_\odot}$'),
          (0.10, 'g7.1.tron', 78, '--', r'$\bm{5\,Z_\odot}$'))
T_JOIN = 4.0                            # log T: F05 below, OPAL above
LOG_NH = (13, 14, 15)                   # n_H [cm^-3] (components plot)
# main plot: a real disk at ~2 ld around 1e7 M_sun (midplane n_H ~ 6e14;
# photosphere, where winds launch, ~1e13 - 1e14)
LOG_NH_REAL = (13, 14, 15)
T_SUB, T_COLOR = 1500.0, 9000.0         # [K] marked
MH = 1.6726e-24
PLOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots',
                    'opacity_vs_T.png')


def read_f05(path=F05):
    """(log T [n], log R [m], log kappa [n, m]) of a Ferguson .tron table"""
    rows, logR = [], None
    for ln in open(path):
        f = ln.split()
        if not f:
            continue
        if f[0] == 'log' and f[1] == 'T':
            logR = np.array(f[2:], float)
        elif logR is not None:
            try:
                rows.append([float(v) for v in f])
            except ValueError:
                pass
    a = np.array(rows)
    o = np.argsort(a[:, 0])
    return a[o, 0], logR, a[o, 1:]


def read_opal(path=OPAL, table=OPAL_TABLE):
    """the same for one OPAL table (short rows: nan past their end)"""
    lines = open(path).read().splitlines()
    k = [i for i, ln in enumerate(lines)
         if ln.startswith(f'TABLE # {table} ') or
         ln.startswith(f'TABLE #{table} ')][-1]
    logR, rows = None, []
    for ln in lines[k + 1:]:
        f = ln.split()
        if not f:
            if rows:
                break
            continue
        if f[0] == 'logT':
            logR = np.array(f[1:], float)
            continue
        if logR is None or f[0] == 'log':
            continue
        v = np.full(len(logR), np.nan)
        vals = [float(x) for x in f[1:]]
        v[:len(vals)] = vals
        rows.append([float(f[0])] + list(v))
    a = np.array(rows)
    return a[:, 0], logR, a[:, 1:]


def kappa_at(tab, logT, logrho):
    """log kappa_R at (log T, log rho) by bilinear interpolation; log R
    clamped to the table (flag: True where clamped)"""
    from scipy.interpolate import RegularGridInterpolator
    lt, lr, lk = tab
    logR = logrho - 3 * (logT - 6)
    clamped = (logR < lr[0]) | (logR > lr[-1])
    f = RegularGridInterpolator((lt, lr), lk, bounds_error=False)
    q = np.column_stack([logT, np.clip(logR, lr[0], lr[-1])])
    return f(q), clamped


# ---- analytic gas opacity components (LTE, Saha) ----
KB_EV, KB_CGS, HPL, CL, ME = 8.6173e-5, 1.3807e-16, 6.6261e-27, 2.9979e10, \
    9.1094e-28
SIG_T = 6.6524e-25
Y_HE = 0.085                    # He / H by number
# electron donors at low T: (abundance / H, ionization energy [eV])
DONORS = ((1.0e-4, 7.8),        # Mg + Si + Fe lumped
          (2.0e-6, 5.14),       # Na
          (1.3e-7, 4.34))       # K
NU = np.logspace(12, 17, 3000)  # [Hz]


def saha_rhs(T, chi, g=1.0):
    """n_{i+1} n_e / n_i [cm^-3] = g (2 pi m_e k T / h^2)^1.5 e^(-chi/kT)"""
    return g * (2 * np.pi * ME * KB_CGS * T / HPL ** 2) ** 1.5 * \
        np.exp(-chi / (KB_EV * T))


def ionization(T, nH, zfac=1.0, y_he=Y_HE):
    """(n_e, n_H0, n_p) [cm^-3] from Saha for H, He (I, II) and the
    DONORS (abundances x zfac; He/H y_he), n_e found by bisection in
    log"""
    def ne_of(ne):
        sH = saha_rhs(T, 13.598)
        xp = sH / (sH + ne)
        s1, s2 = saha_rhs(T, 24.587, 4.0), saha_rhs(T, 54.418, 0.25)
        d = 1 + s1 / ne + s1 * s2 / ne ** 2
        he_e = y_he * (s1 / ne + 2 * s1 * s2 / ne ** 2) / d
        m_e = sum(zfac * a * saha_rhs(T, c) / (saha_rhs(T, c) + ne)
                  for a, c in DONORS)
        return nH * (xp + he_e + m_e), xp
    lo, hi = -30.0, np.log10(nH * 3)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        f, _ = ne_of(10 ** mid)
        lo, hi = (mid, hi) if f > 10 ** mid else (lo, mid)
    ne = 10 ** (0.5 * (lo + hi))
    _, xp = ne_of(ne)
    return ne, nH * (1 - xp), nH * xp


def hminus_bf(lam):
    """H- bound-free cross-section [cm^2] (Wishart 1979 fit, Gray 2005)"""
    a = (1.99654, -1.18267e-5, 2.64243e-6, -4.40524e-10, 3.23992e-14,
         -1.39568e-18, 2.78701e-23)
    s = sum(c * lam ** k for k, c in enumerate(a)) * 1e-18
    return np.where(lam < 16200, np.clip(s, 0, None), 0.0)


def hminus_ff(lam, T):
    """H- free-free per neutral H per unit electron pressure [cm^2 /
    (dyn cm^-2)] (Bell & Berrington 1987 fit, Gray 2005)"""
    x = np.log10(lam)
    f0 = -2.2763 - 1.6850 * x + 0.76661 * x ** 2 - 0.053346 * x ** 3
    f1 = 15.2827 - 9.2846 * x + 1.99381 * x ** 2 - 0.142631 * x ** 3
    f2 = (-197.789 + 190.266 * x - 67.9775 * x ** 2 + 10.6913 * x ** 3
          - 0.625151 * x ** 4)
    lt = np.log10(5040.0 / T)
    return 1e-26 * 10 ** (f0 + f1 * lt + f2 * lt ** 2)


def components(T, nH, zfac=1.0, y_he=Y_HE):
    """monochromatic opacities [cm^-1] on NU of the processes:
    H- bf, H- ff, H bf (n = 1-8, hydrogenic), H ff (Kramers), e scattering"""
    ne, n0, npr = ionization(T, nH, zfac, y_he)
    lam = CL / NU * 1e8                                   # [A]
    stim = -np.expm1(-HPL * NU / (KB_CGS * T))
    th = 5040.0 / T
    pe = ne * KB_CGS * T
    out = {}
    out['hm_bf'] = 4.158e-10 * hminus_bf(lam) * pe * th ** 2.5 * \
        10 ** (0.754 * th) * n0 * stim
    out['hm_ff'] = hminus_ff(lam, T) * pe * n0
    hb = np.zeros_like(NU)
    for n in range(1, 9):
        nu_n = 3.2898e15 / n ** 2
        pop = n0 / 2 * 2 * n ** 2 * np.exp(-13.598 * (1 - 1 / n ** 2)
                                            / (KB_EV * T))
        hb += np.where(NU >= nu_n, 2.815e29 / (n ** 5 * NU ** 3), 0) * pop
    out['h_bf'] = hb * stim
    out['h_ff'] = 3.69e8 * ne * npr * T ** -0.5 * NU ** -3 * stim
    out['es'] = np.full_like(NU, SIG_T * ne)
    return out


def dbdt(T):
    x = HPL * NU / (KB_CGS * T)
    return NU ** 4 * np.exp(-np.clip(x, 0, 700)) / np.clip(
        np.expm1(np.clip(x, 1e-12, 700)), 1e-300, None) ** 2


def rosseland(kappa, T):
    dB = dbdt(T)
    return np.trapz(dB, NU) / np.trapz(dB / np.maximum(kappa, 1e-300), NU)


def rosseland_shares(comp, T):
    """kappa_R of the summed opacity and each process's additive share
    c_i = kappa_R <kappa_i / kappa_tot>, averaged with the Rosseland
    weight dB/dT / kappa_tot (the c_i sum to kappa_R; a process's own
    Rosseland mean is no use: one transparent window sets it to ~0)"""
    tot = np.maximum(sum(comp.values()), 1e-300)
    w = dbdt(T) / tot
    kr = rosseland(tot, T)
    return kr, {k: kr * np.trapz(w * v / tot, NU) / np.trapz(w, NU)
                for k, v in comp.items()}


COMP_STYLE = (('hm_bf', r'$\bm{{\rm H}^-}$\textbf{ bound-free}', 'dodgerblue', '-'),
              ('hm_ff', r'$\bm{{\rm H}^-}$\textbf{ free-free}', 'deepskyblue', ':'),
              ('h_bf', r'\textbf{H bound-free}', 'orangered', '-'),
              ('h_ff', r'\textbf{H free-free}', 'darkorange', ':'),
              ('es', r'\textbf{electron scattering}', '0.5', '--'))
PLOT_COMP = PLOT.replace('.png', '_components.png')


def plot_components(z=0.02):
    """per density (stacked panels): the table kappa_R, the analytic gas
    kappa_R and each process's additive share of it (rosseland_shares;
    gas only: no molecules, dust, metal lines)"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy',
                         'font.size': 24})
    plt.rcParams['text.latex.preamble'] = \
        r'\usepackage{amsmath} \usepackage{bm} \boldmath'
    zz = [m for m in METALS if m[0] == z][0]
    f05 = read_f05(os.path.join(TABLE_DIR, 'f05', 'x', zz[1]))
    opal = read_opal(table=zz[2])
    zfac = z / 0.02
    y_he = (1 - 0.70 - z) / 4 / 0.70          # He/H by number, X = 0.70
    logT = np.linspace(3.0, 5.3, 160)
    T = 10 ** logT
    fig, axs = plt.subplots(len(LOG_NH), 1, sharex=True,
                            figsize=(12, 6.0 * len(LOG_NH)))
    for ax, n in zip(axs, LOG_NH):
        nH = 10.0 ** n
        rho = 1.4 * MH * nH
        lo = logT <= T_JOIN
        k1, _ = kappa_at(f05, logT[lo], np.log10(rho))
        k2, _ = kappa_at(opal, logT[~lo], np.log10(rho))
        ax.loglog(T, 10 ** np.concatenate([k1, k2]), 'k', lw=3,
                  label=r'\textbf{tables (total)}')
        comp = {k: [] for k, *_ in COMP_STYLE}
        tot = []
        for t in T:
            kr, sh = rosseland_shares(components(t, nH, zfac, y_he), t)
            for k in comp:
                comp[k].append(sh[k] / rho)
            tot.append(kr / rho)
        for k, lab, col, ls in COMP_STYLE:
            ax.loglog(T, comp[k], color=col, ls=ls, lw=2.2, label=lab)
        ax.loglog(T, tot, 'k--', lw=1.8, label=r'\textbf{analytic total}')
        for t in (T_SUB, T_COLOR):
            ax.axvline(t, color='0.6', lw=1.2, ls='--')
        ax.text(0.03, 0.92, r'$\bm{n_{\rm H} = 10^{%d}}\,\rm\bf cm^{-3}$' % n
                + r'\textbf{, }' + zz[4],
                transform=ax.transAxes, va='top', fontsize=21)
        ax.set_ylim(1e-5, 3e2)
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', labelsize=19,
                       top=True, right=True)
        ax.tick_params(which='major', length=10, width=2)
        ax.tick_params(which='minor', length=5, width=1.5)
        i = np.argmin(np.abs(T - T_COLOR))
        print(f"log n_H = {n}, T = {T[i]:.0f} K: table {10 ** np.interp(logT[i], logT, np.concatenate([k1, k2])):.3g}; "
              + ', '.join(f"{k} {comp[k][i]:.3g}" for k, *_ in COMP_STYLE)
              + f"; analytic total {tot[i]:.3g} cm^2/g")
    for t, txt in ((T_SUB, r'$\bm{T_{\rm sub}}$'), (T_COLOR,
                                                    r'\textbf{9000\,K}')):
        axs[0].text(t * 1.06, 1.02, txt, fontsize=19, va='bottom',
                    transform=axs[0].get_xaxis_transform())
    axs[len(LOG_NH) // 2].set_ylabel(
        r'$\kappa_{\rm R}$\textbf{ and its shares }$[\rm cm^2\,g^{-1}]$')
    axs[-1].set_xlabel(r'$T\,[\rm K]$')
    axs[-1].set_xlim(T[0], T[-1])
    axs[0].legend(fontsize=15, frameon=False, loc='upper right', ncol=2)
    fig.tight_layout(h_pad=0.2)
    path = PLOT_COMP if z == 0.02 else PLOT_COMP.replace(
        '.png', f'_z{z / 0.02:g}.png')
    fig.savefig(path, dpi=130, bbox_inches='tight')
    print(f"wrote {path}")


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy',
                         'font.size': 24})
    plt.rcParams['text.latex.preamble'] = \
        r'\usepackage{amsmath} \usepackage{bm} \boldmath'
    logT = np.linspace(2.7, 5.3, 400)
    T = 10 ** logT
    fig, ax = plt.subplots(figsize=(12, 8))
    cols = dict(zip(LOG_NH_REAL, ('dodgerblue', 'darkviolet', 'orangered')))
    for z, ffile, otab, ls, zlab in METALS:
        f05 = read_f05(os.path.join(TABLE_DIR, 'f05', 'x', ffile))
        opal = read_opal(table=otab)
        for n in LOG_NH_REAL:
            logrho = np.log10(1.4 * MH) + n
            lo = logT <= T_JOIN
            k1, c1 = kappa_at(f05, logT[lo], logrho)
            k2, c2 = kappa_at(opal, logT[~lo], logrho)
            k, c = np.concatenate([k1, k2]), np.concatenate([c1, c2])
            ax.loglog(T, np.where(c, np.nan, 10 ** k), color=cols[n],
                      lw=2.5, ls=ls,
                      label=(r'$\bm{n_{\rm H} = 10^{%d}}\,\rm\bf cm^{-3}$' % n
                             if z == METALS[0][0] else None))
            ax.loglog(T, np.where(c, 10 ** k, np.nan), color=cols[n], lw=2,
                      ls=':')
            i = np.nanargmax(np.where((T > 5e3) & (T < 3e4), k, -np.inf))
            print(f"Z = {z}, log n_H = {n}: kappa_R max (5e3-3e4 K) "
                  f"{10 ** k[i]:.3g} cm^2/g at {T[i]:.0f} K; 1000 K "
                  f"{10 ** np.interp(3, logT, k):.3g}, 2500 K "
                  f"{10 ** np.interp(np.log10(2500), logT, k):.3g}, 9000 K "
                  f"{10 ** np.interp(np.log10(9000), logT, k):.3g}")
        ax.plot([], [], 'k', ls=ls, lw=2.5, label=zlab)
    for t, txt in ((T_SUB, r'$\bm{T_{\rm sub}}$'), (T_COLOR,
                                                    r'\textbf{9000\,K}')):
        ax.axvline(t, color='0.5', lw=1.2, ls='--')
        ax.text(t * 1.06, 0.97, txt, transform=ax.get_xaxis_transform(),
                va='top', fontsize=21)
    ax.axvline(10 ** T_JOIN, color='0.8', lw=1, ls=':')
    ax.set_xlabel(r'$T\,[\rm K]$')
    ax.set_ylabel(r'$\kappa_{\rm R}\,[\rm cm^2\,g^{-1}]$')
    ax.set_xlim(10 ** logT[0], 10 ** logT[-1])
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', labelsize=19, top=True,
                   right=True)
    ax.tick_params(which='major', length=10, width=2)
    ax.tick_params(which='minor', length=5, width=1.5)
    ax.legend(fontsize=19, frameon=False, loc='lower right')
    fig.tight_layout()
    os.makedirs(os.path.dirname(PLOT), exist_ok=True)
    fig.savefig(PLOT, dpi=130, bbox_inches='tight')
    print(f"wrote {PLOT}")


if __name__ == '__main__':
    import sys
    if 'components' in sys.argv[1:]:
        for m in METALS:
            plot_components(m[0])
    else:
        main()
