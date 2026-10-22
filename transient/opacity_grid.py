#!/usr/bin/env python3
"""opacity_grid.py - Rosseland and Planck mean opacities of gas (and gas +
dust) against temperature from Cloudy, for the disk-atmosphere opacity
question: do the hydrogen-ionization rise near 1e4 K and the dust opacity
below T_sub ~ 1500 K both make opacity-driven structure (FRADO-like puffing,
pillars)?

Each model is one Cloudy zone at a fixed temperature T, density n_H, in
strict thermodynamic equilibrium radiation at T (`blackbody T STE`), so the
gas is in LTE; `save total opacity` gives kappa_nu [cm^-1] (absorption and
scattering, continuum only: Cloudy's lines are not in it, so the Planck mean
is a lower limit where lines matter). Dust: ISM grains whose abundance falls
as exp[-(T_grain/T_sub)^3] (`function sublimation`).

  kappa_R^-1 = int kappa_nu,tot^-1 dB/dT dnu / int dB/dT dnu
  kappa_P    = int kappa_nu,abs B dnu / int B dnu
per unit mass with rho = 1.4 m_H n_H.

Steps (all with --out DIR):
  write   the decks (TEMPS x DENS x {gas, dust})
  run     Cloudy on them (niced, NCPU at a time, unfinished only)
  plot    the means against T (plots/opacity_vs_T.png, opacity_means.npz)

Run:  python3 -m transient.opacity_grid write --out DIR && \
      python3 -m transient.opacity_grid run --out DIR && \
      python3 -m transient.opacity_grid plot --out DIR
Low temperatures with dust (molecules and grains) are the slow models.
"""

import argparse
import glob
import os
import subprocess
import time
from multiprocessing import Pool

import numpy as np

from transient.cloudy_rays import CLOUDY_DIR, EXE, plot_style

TEMPS = np.round(10 ** np.linspace(3, 5, 41))   # [K]
DENS = (10, 12, 14)                             # log n_H [cm^-3]
T_SUB = 1500.0                                  # [K] marked on the plot
T_COLOR = 9000.0                                # [K] Mrk 817 constant comp.
NCPU = 40
MH, KB, HPL, CLIGHT = 1.6726e-24, 1.3807e-16, 6.6261e-27, 2.9979e10
RYD_HZ = 3.28984e15

DECK = '''title opacity T={T:.0f} K, log n_H={n}, {kind}
blackbody t = {T:.0f} K STE
constant temperature {T:.0f} K
hden {n}
{grains}stop zone 1
set dr 0
save total opacity "{name}.opc" last
'''


def model_name(T, n, dust):
    return f"T{T:06.0f}_n{n:02d}_{'dust' if dust else 'gas'}"


def write(out):
    rdir = os.path.join(out, 'models')
    os.makedirs(rdir, exist_ok=True)
    k = 0
    for dust in (False, True):
        for n in DENS:
            for T in TEMPS:
                name = model_name(T, n, dust)
                with open(os.path.join(rdir, name + '.in'), 'w') as fh:
                    fh.write(DECK.format(
                        T=T, n=n, name=name, kind='gas+dust' if dust else
                        'gas', grains='grains ism function sublimation\n'
                        if dust else ''))
                k += 1
    print(f"wrote {k} decks to {rdir}")


def done(rdir, name):
    o, p = (os.path.join(rdir, name + e) for e in ('.out', '.opc'))
    return (os.path.exists(p) and os.path.getsize(p) > 0
            and 'Cloudy exited OK' in open(o, errors='ignore').read()[-500:])


def run_model(args):
    rdir, name = args
    if done(rdir, name):
        return name, 0.0, 'skip'
    t = time.time()
    with open(os.path.join(rdir, name + '.in')) as fi, \
            open(os.path.join(rdir, name + '.out'), 'w') as fo:
        rc = subprocess.run([EXE], stdin=fi, stdout=fo, cwd=rdir,
                            env=dict(os.environ, CLOUDY_DATA_PATH=os.path.join(
                                CLOUDY_DIR, 'data'))).returncode
    return name, time.time() - t, 'ok' if rc == 0 else f'rc={rc}'


def run(out, ncpu=NCPU):
    """Cloudy on every unfinished deck, niced (a long multi-core run);
    fast (high T, gas) models first"""
    os.nice(19)
    rdir = os.path.join(out, 'models')
    names = sorted((os.path.basename(p)[:-3] for p in
                    glob.glob(os.path.join(rdir, '*.in'))),
                   key=lambda n: (n.endswith('dust'), -float(n[1:7])))
    print(f"running {len(names)} models on {ncpu} cores", flush=True)
    with Pool(ncpu) as pool:
        for k, (name, dt, st) in enumerate(pool.imap_unordered(
                run_model, [(rdir, n) for n in names]), 1):
            if st != 'skip':
                print(f"  {k}/{len(names)} {name} {st} {dt:.0f} s",
                      flush=True)


def means(path, T, n):
    """(kappa_R, kappa_P) [cm^2/g] from a .opc file"""
    d = np.loadtxt(path, usecols=(0, 1, 2))
    nu = d[:, 0] * RYD_HZ
    o = np.argsort(nu)
    nu, ktot, kabs = nu[o], d[o, 1], d[o, 2]
    x = HPL * nu / (KB * T)
    B = nu ** 3 / np.expm1(np.clip(x, 1e-300, 700))
    dB = nu ** 4 * np.exp(-np.clip(x, 0, 700)) / np.expm1(
        np.clip(x, 1e-300, 700)) ** 2
    rho = 1.4 * MH * 10 ** n
    kr = np.trapz(dB, nu) / np.trapz(dB / np.maximum(ktot, 1e-300), nu)
    kp = np.trapz(kabs * B, nu) / np.trapz(B, nu)
    return kr / rho, kp / rho


def plot(out):
    rdir = os.path.join(out, 'models')
    res = {}
    for dust in (False, True):
        for n in DENS:
            rows = []
            for T in TEMPS:
                name = model_name(T, n, dust)
                if done(rdir, name):
                    rows.append((T, *means(os.path.join(rdir, name + '.opc'),
                                           T, n)))
            if rows:
                res[(dust, n)] = np.array(rows)
    if not res:
        print("no finished models")
        return
    np.savez(os.path.join(out, 'opacity_means.npz'),
             **{f"{'dust' if d else 'gas'}_n{n}": v
                for (d, n), v in res.items()})
    plt = plot_style()
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 12))
    cols = dict(zip(DENS, ('dodgerblue', 'darkviolet', 'orangered')))
    for ax, j, lab in ((axs[0], 1, r'$\kappa_{\rm R}$'),
                       (axs[1], 2, r'$\kappa_{\rm P}$')):
        for (dust, n), v in sorted(res.items()):
            ax.loglog(v[:, 0], v[:, j], color=cols[n], lw=2,
                      ls='--' if dust else '-',
                      label=(r'$\bm{\log n_{\rm H} = %d}$' % n
                             + (r'\textbf{, dust}' if dust else '')))
        for t, txt in ((T_SUB, r'$\bm{T_{\rm sub}}$'),
                       (T_COLOR, r'\textbf{9000\,K}')):
            ax.axvline(t, color='0.5', lw=1, ls=':')
        ax.set_ylabel(lab + r'$\,[\rm cm^2\,g^{-1}]$')
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', labelsize=19,
                       top=True, right=True)
        ax.tick_params(which='major', length=10, width=2)
        ax.tick_params(which='minor', length=5, width=1.5)
    axs[0].text(T_SUB * 1.05, 0.95, r'$\bm{T_{\rm sub}}$', fontsize=19,
                transform=axs[0].get_xaxis_transform(), va='top')
    axs[0].text(T_COLOR * 1.05, 0.95, r'\textbf{9000\,K}', fontsize=19,
                transform=axs[0].get_xaxis_transform(), va='top')
    axs[0].legend(fontsize=15, frameon=False, ncol=2, loc='lower right')
    axs[-1].set_xlabel(r'$T\,[\rm K]$')
    fig.tight_layout(h_pad=0.2)
    os.makedirs(os.path.join(out, 'plots'), exist_ok=True)
    path = os.path.join(out, 'plots', 'opacity_vs_T.png')
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {path}: " + ', '.join(f"{'dust' if d else 'gas'} n{n} "
                                        f"{len(v)} T" for (d, n), v in
                                        sorted(res.items())))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('step', choices=['write', 'run', 'plot', 'all'])
    ap.add_argument('--out', default='/data2/jhuang/runs/opacity_grid')
    ap.add_argument('--ncpu', type=int, default=NCPU)
    a = ap.parse_args()
    if a.step in ('write', 'all'):
        write(a.out)
    if a.step in ('run', 'all'):
        run(a.out, a.ncpu)
    if a.step in ('plot', 'all'):
        plot(a.out)


if __name__ == '__main__':
    main()
