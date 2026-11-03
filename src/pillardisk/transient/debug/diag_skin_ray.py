#!/usr/bin/env python3
"""diag_skin_ray.py - where does a Cloudy ray's reflected ionising light come
from: the hot thin gas the ray crosses before the dense skin, or the skin?
(docs/NOTE_diffuse_field.md)

Each test ray (TEST_RAYS, from train/h0.0) is rerun with its deck starting
at the skin (the first deck point with n_H >= N_SURF), phi(H) diluted to
there (lphi_skin), everything else unchanged. If the reflected ionising
fraction stays ~1, it comes from the skin (recombination photons escaping
its lit face); if it drops, from the hot gas in front.

Steps:
  write   decks <OUT>/<ray>_skin.in
  run     Cloudy on them (one core each, nice 19; ~10-40 min per ray)
  report  reflected / incident ionising photons and energy per band, the
          original rays against the skin rays: printed, saved to
          <OUT>/skin_test.npz, plotted to <OUT>/skin_test.png
          (reflected / incident nu F_nu against wavelength)
Run:  python3 -m pillardisk.transient.debug.diag_skin_ray all
"""
import argparse
import os
import shutil
import sys

import numpy as np

from pillardisk.transient.cloudy_rays import N_SURF, run_one
from pillardisk.transient.ray_emulator import read_deck, ray_features
from pillardisk.transient.sim_spectrum import CONFIG, LD_CM

SRC = '/data2/jhuang/runs/cloudy_rays/train/h0.0/rays'
OUT = '/data2/jhuang/runs/cloudy_rays/diffuse/skin_test'
# two rays entering through thin hot gas (n_H 2e9 and 6e8 at the entry,
# Te 1e6 and 3e6 K; reflected ionising energy 0.96, 0.94). Rays entering
# dense gas reflect as much with nothing to cut: ray_p215_a06 (n_H 1.3e10,
# Te 2e5 K at its face) 0.94, ray_p044_a07 (3e11, 4e4 K) 0.70
TEST_RAYS = ('ray_p178_a09', 'ray_p009_a00')
H, C = 6.626e-27, 2.998e10
BANDS = ((504.3, 911.76, 'H 13.6-24.6 eV'), (227.8, 504.3, 'He I 24.6-54.4'),
         (50.0, 227.8, 'He II 54-250 eV'), (0.0, 50.0, '> 250 eV'))


def skin_deck(name):
    """the deck of ray name cut to start at its skin"""
    lphi, s, n = read_deck(os.path.join(SRC, name + '.in'))
    f, sk = ray_features(lphi, s, n, CONFIG['R0_LD'] * LD_CM)
    lsk = np.log10(sk)
    lines = open(os.path.join(SRC, name + '.in')).read().split('\n')
    out, on, pts = [], False, []
    for ln in lines:
        if ln.startswith('dlaw table'):
            on = True
            out.append(ln)
            continue
        if on and ln.startswith('end of dlaw'):
            # the original's padding: a point 0.01 dex before the first
            # real one (added here) and one 0.01 dex past the last (kept)
            real = [p for p in pts if p[0] >= lsk - 1e-9]
            body = [(real[0][0] - 0.01, real[0][1])] + real
            out += [f'continue {r:.7f} {v:.4f}' for r, v in body]
            out.append(ln)
            on = False
            continue
        if on and ln.startswith('continue'):
            pts.append(tuple(float(v) for v in ln.split()[1:3]))
            continue
        if ln.startswith('phi(H)'):
            ln = f"phi(H) {f['lphi_skin']:.4f}"
        elif ln.startswith('radius'):
            ln = f'radius {lsk:.5f} {ln.split()[2]}'
        elif ln.startswith('title'):
            ln += ' - cut to start at the skin'
        out.append(ln.replace(f'"{name}.', f'"{name}_skin.'))
    return '\n'.join(out), f['lphi_skin'], lsk


def write():
    os.makedirs(OUT, exist_ok=True)
    for name in TEST_RAYS:
        deck, lphi, lsk = skin_deck(name)
        path = os.path.join(OUT, name + '_skin.in')
        with open(path, 'w') as fh:
            fh.write(deck)
        _, s, n = read_deck(path)
        print(f'{path}: starts at log r {lsk:.4f}, phi(H) {lphi:.3f}, '
              f'n_H {n[0]:.3g}..{n.max():.3g}, {len(s)} points')
        # the originals' outputs, for the report
        for ext in ('in', 'conA'):
            shutil.copy(os.path.join(SRC, f'{name}.{ext}'),
                        os.path.join(OUT, f'{name}.{ext}'))


def run():
    from multiprocessing import Pool
    os.nice(19)
    with Pool(len(TEST_RAYS)) as p:
        for name, dt, st in p.imap_unordered(
                run_one, [(OUT, n + '_skin') for n in TEST_RAYS]):
            print(f'{name}: {st} in {dt / 60:.1f} min', flush=True)


def bands(path):
    """per BANDS: (incident, reflected) photon rates and energies, and the
    spectrum (wave, incident, reflected cont + lines)"""
    import pandas as pd
    a = pd.read_csv(path, sep='\t', usecols=[0, 1, 5, 7], skiprows=1,
                    header=None).values[::-1]
    dl = np.gradient(np.log(a[:, 0]))
    ph = a[:, 0] * 1e-8 / (H * C) * dl
    refl = a[:, 2] + a[:, 3]
    rows = []
    for lo, hi, _ in BANDS:
        m = (a[:, 0] >= lo) & (a[:, 0] < hi)
        rows.append((np.sum(a[m, 1] * ph[m]), np.sum(refl[m] * ph[m]),
                     np.sum(a[m, 1] * dl[m]), np.sum(refl[m] * dl[m])))
    return np.array(rows), (a[:, 0], a[:, 1], refl)


def report():
    from pillardisk.transient.cloudy_rays import plot_style
    res, spec = {}, {}
    for name in TEST_RAYS:
        for tag in ('', '_skin'):
            p = os.path.join(OUT, f'{name}{tag}.conA')
            if not os.path.exists(p) or os.path.getsize(p) == 0:
                print(f'{p}: missing (run first)')
                continue
            res[name + tag], spec[name + tag] = bands(p)
    for name in TEST_RAYS:
        print(f'\n{name}: reflected / incident (original -> skin start)')
        print(f"{'band':18s} {'photons':>16s} {'energy':>16s}")
        for b, (_, _, lab) in enumerate(BANDS):
            cells = []
            for k in (1, 3):
                v = [res[t][b, k] / res[t][b, k - 1] if t in res else np.nan
                     for t in (name, name + '_skin')]
                cells.append(f'{v[0]:6.2f} -> {v[1]:5.2f}')
            print(f'{lab:18s} {cells[0]:>16s} {cells[1]:>16s}')
        v = [res[t][:, 1].sum() / res[t][:, 0].sum() if t in res else np.nan
             for t in (name, name + '_skin')]
        print(f"{'all ionising':18s} {v[0]:6.2f} -> {v[1]:5.2f} (photons)")
    np.savez(os.path.join(OUT, 'skin_test.npz'),
             bands=np.array([b[2] for b in BANDS]),
             **{f'bands_{k}': v for k, v in res.items()},
             **{f'spec_{k}': np.array(v) for k, v in spec.items()})
    plt = plot_style()
    fig, axs = plt.subplots(len(TEST_RAYS), 1, figsize=(11, 5.2 *
                                                         len(TEST_RAYS)),
                            sharex=True, sharey=True)
    for ax, name in zip(np.atleast_1d(axs), TEST_RAYS):
        for tag, col, lab in (('', 'k', 'original'),
                              ('_skin', 'crimson', 'skin start')):
            if name + tag in spec:
                w, inc, refl = spec[name + tag]
                m = (w < 912) & (inc > 0)
                ax.plot(w[m], refl[m] / inc[m], color=col, lw=1.2,
                        label=rf'\textbf{{{lab}}}')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-3, 30)
        ax.set_ylabel(r'\textbf{reflected / incident}')
        ax.text(0.03, 0.93, rf"\textbf{{{name.replace('_', ' ')}}}",
                transform=ax.transAxes, va='top', fontsize=18)
        ax.tick_params(which='both', direction='in', labelsize=17)
        ax.legend(fontsize=16, loc='lower right')
    np.atleast_1d(axs)[-1].set_xlabel(r'$\bm{\lambda\,[\rm\AA]}$')
    fig.tight_layout()
    path = os.path.join(OUT, 'skin_test.png')
    fig.savefig(path, dpi=120)
    print(f'\nwrote {path} and skin_test.npz')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('step', choices=['write', 'run', 'report', 'all'])
    a = ap.parse_args()
    if a.step in ('write', 'all'):
        write()
    if a.step in ('run', 'all'):
        run()
    if a.step in ('report', 'all'):
        report()


if __name__ == '__main__':
    sys.exit(main())
