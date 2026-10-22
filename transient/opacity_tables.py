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
"""

import os

import numpy as np

TABLE_DIR = '/data2/jhuang/data/opacity_tables'
F05 = os.path.join(TABLE_DIR, 'f05', 'x', 'g7.02.tron')
OPAL = os.path.join(TABLE_DIR, 'opal', 'GN93hz')
OPAL_TABLE = 73                         # X = 0.70, Z = 0.02
T_JOIN = 4.0                            # log T: F05 below, OPAL above
LOG_NH = (10, 12, 14)                   # n_H [cm^-3]
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


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy',
                         'font.size': 24})
    plt.rcParams['text.latex.preamble'] = \
        r'\usepackage{amsmath} \usepackage{bm} \boldmath'
    f05, opal = read_f05(), read_opal()
    logT = np.linspace(f05[0][0], 5.3, 400)
    fig, ax = plt.subplots(figsize=(12, 8))
    cols = dict(zip(LOG_NH, ('dodgerblue', 'darkviolet', 'orangered')))
    for n in LOG_NH:
        logrho = np.log10(1.4 * MH) + n
        lo = logT <= T_JOIN
        k1, c1 = kappa_at(f05, logT[lo], logrho)
        k2, c2 = kappa_at(opal, logT[~lo], logrho)
        k, c = np.concatenate([k1, k2]), np.concatenate([c1, c2])
        T = 10 ** logT
        ax.loglog(T, np.where(c, np.nan, 10 ** k), color=cols[n], lw=2.5,
                  label=r'$\bm{n_{\rm H} = 10^{%d}}\,\rm\bf cm^{-3}$' % n)
        ax.loglog(T, np.where(c, 10 ** k, np.nan), color=cols[n], lw=2,
                  ls=':')
        i = np.nanargmax(np.where((T > 5e3) & (T < 3e4), k, -np.inf))
        print(f"log n_H = {n}: kappa_R peak near the H-ionization bump "
              f"{10 ** k[i]:.3g} cm^2/g at T = {T[i]:.0f} K; at 1000 K "
              f"{10 ** np.interp(3, logT, k):.3g}, at 2500 K "
              f"{10 ** np.interp(np.log10(2500), logT, k):.3g}")
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
    main()
