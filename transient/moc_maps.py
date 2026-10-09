#!/usr/bin/env python3
"""moc_maps.py - maps of a finished MOCASSIN run of the disk snapshot:
ion fractions and temperature (from the grid files) and line emission
(from mocassinPlot) in three views.

Reads, from <run>/output:
  grid0.out  axes [cm] and per-cell (active index, converged, black)
  grid1.out  per-cell Te [K], Ne, n_H [cm^-3]
  grid2.out  per-cell ion fractions, one line per element switched on
  plot.out   (optional) per-cell line luminosities [1e36 erg/s] from
             mocassinPlot with the line list in <run>/input/plot.in; the
             columns are named by <run>/input/plot_lines.txt (written with
             plot.in by transient.moc_plotin), else they must be the
             original seven: H alpha, H beta, C IV 1551+1548,
             Mg II 2804+2796, [O III] 5007
and <run>/grid.npz (written by transient.moc_convert) for the r0 scale.

Views (one figure each):
  faceon    face-on projection: line surface brightness [erg/s/cm^2]
            from the observer side only (z > 0, half the midplane cell); ion
            fractions and Te n_H-weighted along z; column density N_H
  midplane  z = 0 slice
  vertical  y = 0 slice (x-z plane through the lamp)
Panels: Te, x(H+), x(C3+), x(Mg+) / n_H (N_H face-on), H alpha, C IV,
Mg II.

Also prints a consistency check of the summed plot.out lines against
the volume-integrated totals in output/lineFlux.out.

Outputs: <run>/plots/moc_maps_{faceon,midplane,vertical}.png
Run:  python3 -m transient.moc_maps [--run /data2/jhuang/runs/mocassin/moc_coarse]
"""

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

DEFAULT_RUN = '/data2/jhuang/runs/mocassin/moc_coarse'
NSTAGES = 7
# plot.out columns of the original seven-line plot.in (used when there is
# no input/plot_lines.txt); the first seven of moc_plotin's list too
PLOT_COLS = ['Ha', 'Hb', 'CIV1551', 'CIV1548', 'MgII2804', 'MgII2796',
             'OIII5007']
PLOT_CODES = [1, 2, 1393, 1394, 6884, 6885, 3741]


def read_line_table(run):
    """input/plot_lines.txt as a list of dicts (None if absent)."""
    path = os.path.join(run, 'input', 'plot_lines.txt')
    if not os.path.exists(path):
        return None
    rows = []
    for line in open(path):
        if line.startswith('#'):
            continue
        t = line.split()
        rows.append(dict(name=t[0], code=int(t[1]), lam=float(t[2]),
                         air=int(t[3]), mass=float(t[4]),
                         ratio=float(t[5])))
    return rows


def read_tokens(path):
    with open(path) as fh:
        return np.array(fh.read().split(), dtype=float)


def load_run(run):
    out = os.path.join(run, 'output')
    with open(os.path.join(out, 'grid0.out')) as fh:
        tok = fh.read().split()
    nx, ny, nz = (int(t) for t in tok[1:4])
    p = 7
    ax = []
    for n in (nx, ny, nz):
        ax.append(np.array(tok[p:p + n], dtype=float))
        p += n
    cell = np.array(tok[p:p + 3 * nx * ny * nz], dtype=float)
    cell = cell.reshape(nx, ny, nz, 3)        # z fastest, as written
    shape = (nx, ny, nz)
    g1 = read_tokens(os.path.join(out, 'grid1.out')).reshape(*shape, 3)

    # elements switched on = nonzero abundance; per element min(Z+1, 7)
    # stages, in order of Z
    abun = [float(line.split()[0]) for line in
            open(os.path.join(run, 'input', 'abun.in'))]
    on = [z for z in range(1, len(abun) + 1) if abun[z - 1] > 0]
    nst = [min(z + 1, NSTAGES) for z in on]
    g2 = read_tokens(os.path.join(out, 'grid2.out'))
    g2 = g2.reshape(*shape, sum(nst))
    ions, k = {}, 0
    for z, n in zip(on, nst):
        ions[z] = g2[..., k:k + n]
        k += n

    # "black" cells failed the T or ionization balance in the last update
    # (update_mod.f90); outputGas leaves them out of the line totals, so
    # mask them with the empty cells
    d = dict(axes=ax, active=(cell[..., 0] > 0) & (cell[..., 2] == 0),
             conv=cell[..., 1] > 0, black=cell[..., 2] != 0,
             Te=g1[..., 0], Ne=g1[..., 1], nH=g1[..., 2], ions=ions)

    pfile = os.path.join(out, 'plot.out')
    if os.path.exists(pfile):
        pl = read_tokens(pfile).reshape(nx * ny * nz, -1)[:, 1:]
        table = read_line_table(run)
        if table is not None and len(table) == pl.shape[1]:
            names = [r['name'] for r in table]
            d['line_table'] = table
        elif pl.shape[1] == len(PLOT_COLS):
            names = PLOT_COLS
            d['line_table'] = [dict(name=n, code=c) for n, c in
                               zip(PLOT_COLS, PLOT_CODES)]
        else:
            raise SystemExit(f'plot.out has {pl.shape[1]} lines, matching '
                             'neither input/plot_lines.txt nor the original '
                             'seven: rerun mocassinPlot')
        d['lines'] = {name: pl[:, i].reshape(shape)
                      for i, name in enumerate(names)}

    # r0 in cm: grid.npz stores the same cell centres in units of r0
    npz = np.load(os.path.join(run, 'grid.npz'))
    i = np.argmax(np.abs(npz['x']))
    d['r0_cm'] = ax[0][i] / npz['x'][i]
    # gas cells MOCASSIN emptied because they lie inside Rin
    d['rin_cut'] = (npz['nH'] > 0) & (cell[..., 0] <= 0)
    return d


def check_totals(d, run):
    """Summed plot.out luminosities vs the lineFlux.out totals (the
    Formal Solution block lists every line as '... ratio code'): checks
    that each plot.out column is the line its code claims."""
    txt = open(os.path.join(run, 'output', 'lineFlux.out')).read()
    hb_ref = float(re.search(r'Hbeta \[E36 erg/s\]:\s+(\S+)', txt).group(1))
    hb = d['lines']['Hb'].sum()
    print(f"  H beta: plot.out {hb:.5g}, lineFlux.out {hb_ref:.5g} "
          f"[1e36 erg/s]")
    fs = txt[txt.index('Formal Solution'):txt.index('Component:')]
    dev = []
    for r in d['line_table']:
        m = re.search(rf'\s(\S+)\s+{r["code"]}\s*\n', fs)
        ref = float(m.group(1)) if m else np.nan
        val = d['lines'][r['name']].sum() / hb
        dev.append((abs(val / ref - 1) if ref > 0 else np.inf, r['name'],
                    val, ref))
    dev.sort(reverse=True)
    print(f"  {len(dev)} lines vs lineFlux.out (ratio to H beta): max "
          f"deviation {100 * dev[0][0]:.2f}% ({dev[0][1]}: {dev[0][2]:.4g} "
          f"vs {dev[0][3]:.4g}); median {100 * np.median([x[0] for x in dev]):.2f}%")


def edges(c):
    m = 0.5 * (c[1:] + c[:-1])
    return np.concatenate([[2 * c[0] - m[0]], m, [2 * c[-1] - m[-1]]])


def view_maps(d, view):
    """2D maps for one view; returns dict of arrays and the axes."""
    act, nH = d['active'], d['nH']
    w = np.where(act, nH, 0.0)
    ions = d['ions']
    qty = {'xH': ions[1][..., 1], 'xC4': ions[6][..., 3],
           'xMg2': ions[12][..., 1], 'Te': d['Te']}
    lines = {}
    if 'lines' in d:
        L = d['lines']
        lines = {'Ha': L['Ha'], 'CIV': L['CIV1551'] + L['CIV1548'],
                 'MgII': L['MgII2804'] + L['MgII2796']}
    x, y, z = (a / d['r0_cm'] for a in d['axes'])
    out = {}
    if view == 'faceon':
        wsum = w.sum(axis=2)
        zc = d['axes'][2]
        zside = np.where(zc > 0, 1.0, 0.0)
        zside[np.argmin(np.abs(zc))] = 0.5
        area = np.outer(np.diff(edges(d['axes'][0])),
                        np.diff(edges(d['axes'][1])))       # [cm^2]
        for k, v in qty.items():
            out[k] = np.where(wsum > 0, (v * w).sum(axis=2)
                              / np.where(wsum > 0, wsum, 1), np.nan)
        for k, v in lines.items():
            # observer side only (z > 0, plus half the midplane cell): the
            # far side is hidden behind the N_H ~ 1e24-25 cm^-2 midplane,
            # as in the pipeline maps; then per unit projected area
            out[k] = (v * zside).sum(axis=2) * 1e36 / area   # erg/s/cm^2
        dz = np.diff(edges(d['axes'][2]))               # [cm]
        out['nH'] = (w * dz).sum(axis=2)                # column N_H [cm^-2]
        # columns that lost gas to Rin (the staircase around the inner
        # hole) miss their densest cells: blank them, not show them dim
        cut = d['rin_cut'].any(axis=2)
        for k in out:
            out[k] = np.where(cut, np.nan, out[k])
        h, v_ = x, y
        labels = (r'$x\,[r_0]$', r'$y\,[r_0]$')
    else:
        if view == 'midplane':
            sl = (slice(None), slice(None), len(z) // 2)
            h, v_ = x, y
            labels = (r'$x\,[r_0]$', r'$y\,[r_0]$')
        else:
            sl = (slice(None), len(y) // 2, slice(None))
            h, v_ = x, z
            labels = (r'$x\,[r_0]$', r'$z\,[r_0]$')
        a = act[sl]
        for k, v in qty.items():
            out[k] = np.where(a, v[sl], np.nan)
        for k, v in lines.items():
            out[k] = v[sl]
        out['nH'] = np.where(a, nH[sl], np.nan)
    return out, h, v_, labels


def add_colorbar(pc, ax):
    """Colorbar with exactly the height of the parent axes (same as
    transient.sim_spectrum.add_colorbar, which cannot be imported on the
    server: it needs athena_read)."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes('right', size='4.5%',
                                              pad=0.12)
    cb = ax.figure.colorbar(pc, cax=cax)
    cb.ax.minorticks_on()
    cb.ax.tick_params(which='major', direction='in', length=8, width=1.5,
                      labelsize=13)
    cb.ax.tick_params(which='minor', direction='in', length=4, width=1.0)
    return cb


def panels(view):
    """Row 1: Te and ion fractions; row 2: density and line emission."""
    dens = (r'$N_{\rm H}\,[\rm cm^{-2}]$' if view == 'faceon'
            else r'$n_{\rm H}\,[\rm cm^{-3}]$')
    # face-on: observer-side surface brightness; slices: per-cell L
    line = (r'$S(\rm %s)\,[\rm erg\,s^{-1}\,cm^{-2}]$' if view == 'faceon'
            else r'$L(\rm %s)\,[10^{36}\,\rm erg\,s^{-1}]$')
    return [('Te', r'$T_e\,[\rm K]$', 'logT'),
            ('xH', r'$x(\rm H^+)$', 'lin01'),
            ('xC4', r'$x(\rm C^{3+})$', 'lin01'),
            ('xMg2', r'$x(\rm Mg^{+})$', 'lin01'),
            ('nH', dens, 'logn'),
            ('Ha', line % r'H\alpha', 'loglum'),
            ('CIV', line % r'C\,IV\,1549', 'loglum'),
            ('MgII', line % r'Mg\,II\,2798', 'loglum')]


def plot_view(maps, h, v, labels, title, path, view):
    he, ve = edges(h), edges(v)
    aspect = (ve[-1] - ve[0]) / (he[-1] - he[0])
    # map width per panel ~4.4 in (24 in / 4 minus colorbar and labels);
    # rows get that times the aspect plus ~1.4 in for title and x label
    fig, axs = plt.subplots(2, 4, figsize=(24, 2 * (4.4 * aspect + 1.4) + 0.6))
    # one shared colour range for all line maps, so they compare directly
    lmax = max([np.nanmax(maps[k]) for k, _, kind in panels(view)
                if kind == 'loglum' and k in maps] or [1.0])
    lnorm = LogNorm(lmax * 1e-5, lmax)
    for n, (ax, (key, lab, kind)) in enumerate(zip(axs.ravel(), panels(view))):
        if key not in maps:
            ax.text(0.5, 0.5, r'$\rm run\ mocassinPlot$', ha='center',
                    va='center', transform=ax.transAxes)
            ax.set_axis_off()
            continue
        m = maps[key].T
        if kind == 'lin01':
            im = ax.pcolormesh(he, ve, m, vmin=0, vmax=1, cmap='viridis')
        elif kind == 'logT':
            im = ax.pcolormesh(he, ve, m, norm=LogNorm(5e2, 1e5),
                               cmap='inferno')
        elif kind == 'logn':
            m = np.where(m > 0, m, np.nan)
            im = ax.pcolormesh(he, ve, m, norm=LogNorm(np.nanmin(m),
                                                       np.nanmax(m)),
                               cmap='viridis')
        else:
            m = np.where(m > 0, m, np.nan)
            im = ax.pcolormesh(he, ve, m, norm=lnorm, cmap='magma')
        ax.set_aspect('equal')
        ax.set_title(lab)
        ax.set_xlabel(labels[0])
        if n % 4 == 0:                       # y label on the left column only
            ax.set_ylabel(labels[1])
        add_colorbar(im, ax)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default=DEFAULT_RUN)
    a = ap.parse_args()
    d = load_run(a.run)
    nact = d['active'].sum()
    print(f"{a.run}: grid {d['active'].shape}, {nact} active cells, "
          f"{(d['conv'] & d['active']).sum() / nact:.1%} converged, "
          f"r0 = {d['r0_cm']:.4g} cm")
    hot = d['active'] & (d['Te'] > 2.6e4)
    print(f"  cells with Te > 26000 K: {hot.sum()} "
          f"(n_H median {np.median(d['nH'][hot]) if hot.any() else 0:.3g})")
    if 'lines' in d:
        check_totals(d, a.run)
    else:
        print('  no output/plot.out yet: line panels left empty')
    pdir = os.path.join(a.run, 'plots')
    os.makedirs(pdir, exist_ok=True)
    titles = {'faceon': r'$\rm face\mbox{-}on\ (lines:\ observer\ side\ z>0;\ x,\ T_e:\ n_H\mbox{-}weighted)$',
              'midplane': r'$\rm midplane\ slice\ (z=0)$',
              'vertical': r'$\rm vertical\ slice\ (y=0)$'}
    for view in ('faceon', 'midplane', 'vertical'):
        maps, h, v, labels = view_maps(d, view)
        path = os.path.join(pdir, f'moc_maps_{view}.png')
        plot_view(maps, h, v, labels, titles[view], path, view)
        print(f'  wrote {path}')


if __name__ == '__main__':
    main()
