#!/usr/bin/env python3
"""ray_emulator.py - emulate the per-ray Cloudy outputs of a cloudy_rays run
from each ray's inputs, so that finer ray sets can be made without Cloudy.

Inputs (features, from the ray's deck alone: what any new ray traced through
the sim would have): the incident flux at the entry, the sim's n_H along the
ray (dlaw table). The ray "lands" where n_H first exceeds N_SURF:
  lphi_skin  log flux there (geometric 1/s^2 dilution from the entry)
  lN_pre     log column [cm^-2] from the entry to there
  lL_skin    log density scale length there [r0] (n / dn/ds)
  ln_p1, ln_p3   log n_H 0.01 and 0.03 r0 past it
  lN_post    log column from there to the ray's end
  lnmax      log of the highest n_H along the ray
  sin_graze  sine of the angle between the ray and the skin there (|ray .
             grad n_H| / |grad n_H|, from the sim): unlike the direction
             from the lamp, comparable between lamp heights
  R_skin, z_skin   the landing point [r0]
  lands      1 if n_H ever exceeds N_SURF (else the skin is the n_H peak)
  h          the lamp height [r0]
Targets (from the run's ray cache), per line (Ha, Hb, C IV, Mg II, He II):
  lL_<line>  log sum 4 pi j s^2 ds [erg/s/sr x 4 pi]
  dc_<line>  its emissivity-weighted point minus the skin [r0]
  lw_<line>  log of its layer (5-95% of the emission in depth) [r0]

`validate` holds out whole phi columns (GROUP_FOLDS folds; neighbouring rays
in a column are correlated) and compares, ray by ray and per held-out column
total, the emulator (one gradient-boosted tree model per target) with a
baseline that interpolates the two neighbouring training columns in alpha
and phi. `heights` trains on all but one lamp height (HEIGHT_RUNS) and
predicts the one left out: can a height be emulated from the others?

Run:  python3 -m transient.ray_emulator validate --out RUN_DIR
      python3 -m transient.ray_emulator heights
(the tree models use EMU_THREADS threads: all cores fight on a busy machine)
"""

import argparse
import os

EMU_THREADS = 8
os.environ.setdefault('OMP_NUM_THREADS', str(EMU_THREADS))

import numpy as np                                        # noqa: E402

from transient.cloudy_rays import (load_rays, N_SURF, plot_style,  # noqa
                                   density_interp, DUMP)
from transient.sim_spectrum import (CONFIG, LD_CM, load_sim,  # noqa: E402
                                    to_physical)

LINES = dict(ha=(1,), hb=(2,), c4=(3, 4), mg2=(5, 6), he2=(7,))
LINE_TEX = dict(ha=r'H$\bm{\alpha}$', hb=r'H$\bm{\beta}$', c4=r'C\,IV',
                mg2=r'Mg\,II', he2=r'He\,II\,1640')
FEATURES = ('lphi_skin', 'lN_pre', 'lL_skin', 'ln_p1', 'ln_p3', 'lN_post',
            'lnmax', 'sin_graze', 'R_skin', 'z_skin', 'lands', 'h')
RUNS = '/data2/jhuang/runs/cloudy_rays'
HEIGHT_RUNS = (f'{RUNS}/train/h0.0', f'{RUNS}/train/h0.2',
               f'{RUNS}/edd10_h0.5_adapt', f'{RUNS}/train/h1.0')
GROUP_FOLDS = 8


def read_deck(path):
    """(log phi(H) at the entry, s [cm], n_H [cm^-3]) from a ray's deck,
    the dlaw padding points dropped"""
    lphi, pts, on = None, [], False
    for ln in open(path):
        if ln.startswith('phi(H)'):
            lphi = float(ln.split()[1])
        elif ln.startswith('dlaw table'):
            on = True
        elif ln.startswith('end of dlaw'):
            on = False
        elif on and ln.startswith('continue'):
            pts.append([float(v) for v in ln.split()[1:3]])
    a = np.array(pts)[1:-1]
    return lphi, 10 ** a[:, 0], 10 ** a[:, 1]


def ray_features(lphi, s, n, alpha, h, r0):
    """the FEATURES of one ray (see the module docstring)"""
    col = np.concatenate([[0], np.cumsum(0.5 * (n[1:] + n[:-1])
                                         * np.diff(s))])
    hit = np.nonzero(n >= N_SURF)[0]
    k = hit[0] if len(hit) else int(np.argmax(n))
    sk = s[k]
    lnn = np.log(n)
    j0, j1 = max(k - 3, 0), min(k + 3, len(s) - 1)
    g = (lnn[j1] - lnn[j0]) / max(s[j1] - s[j0], 1.0)  # d ln n / ds [1/cm]
    lsc = np.log10(np.clip(1 / max(abs(g), 1e-30) / r0, 1e-6, 10))
    at = lambda ds: np.log10(np.interp(sk + ds * r0, s, n))
    return dict(lphi_skin=lphi + 2 * np.log10(s[0] / sk),
                lN_pre=np.log10(col[k] + 1e10),
                lL_skin=lsc, ln_p1=at(0.01), ln_p3=at(0.03),
                lN_post=np.log10(col[-1] - col[k] + 1e10),
                lnmax=np.log10(n.max()), alpha=alpha,
                R_skin=sk * np.sin(alpha) / r0,
                z_skin=h + sk * np.cos(alpha) / r0,
                lands=float(len(hit) > 0)), sk


def sin_graze(sim, itp, phi, alpha, h, sk_r0, eps=0.005):
    """|ray . grad n_H| / |grad n_H| at the landing point (sk_r0 [r0]
    from the lamp at z = h along (alpha, phi)), log n_H differenced over
    +- eps r0 in x, y, z"""
    d = np.array([np.sin(alpha) * np.cos(phi), np.sin(alpha) * np.sin(phi),
                  np.cos(alpha)])
    p = np.array([0.0, 0.0, h]) + sk_r0 * d
    pts = p[None, :] + eps * np.concatenate([np.eye(3), -np.eye(3)])
    r = np.linalg.norm(pts, axis=1)
    q = np.column_stack([np.mod(np.arctan2(pts[:, 1], pts[:, 0]), 2 * np.pi),
                         np.clip(np.arccos(pts[:, 2] / r), sim['th'][0],
                                 sim['th'][-1]),
                         np.clip(r, sim['r'][0], sim['r'][-1])])
    v = itp(q)
    g = (v[:3] - v[3:]) / (2 * eps)
    gn = np.linalg.norm(g)
    return float(abs(g @ d) / gn) if gn > 0 else 0.0


def ray_targets(e, s0, sk, r0):
    """the targets of one ray from its .ems (depth, 4 pi j columns)"""
    d = e[:, 0]
    s = s0 + d
    ds = np.gradient(d)
    out = {}
    for key, cols in LINES.items():
        w = e[:, list(cols)].sum(axis=1) * s ** 2 * ds
        tot = w.sum()
        if not tot > 0:
            out.update({f'lL_{key}': np.nan, f'dc_{key}': np.nan,
                        f'lw_{key}': np.nan})
            continue
        cw = np.cumsum(w) / tot
        q5, q95 = np.interp([0.05, 0.95], cw, s)
        out[f'lL_{key}'] = np.log10(tot)
        out[f'dc_{key}'] = (np.sum(w * s) / tot - sk) / r0
        out[f'lw_{key}'] = np.log10(max(q95 - q5, 1e-6 * r0) / r0)
    return out


def build_table(out):
    """features and targets of every finished ray of the run (adaptive
    layout: rays_meta.npz), cached as <out>/emulator/table.npz"""
    path = os.path.join(out, 'emulator', 'table.npz')
    if os.path.exists(path):
        d = np.load(path)
        if all(k in d.files for k in FEATURES):
            return {k: d[k] for k in d.files}
    sim = load_sim(DUMP)
    itp = density_interp(sim, to_physical(sim, CONFIG))
    m = np.load(os.path.join(out, 'rays_meta.npz'))
    h, r0 = float(m['hlamp_r0']), CONFIG['R0_LD'] * LD_CM
    meta = {str(n): (int(p), float(a), float(d0)) for n, p, a, d0 in
            zip(m['name'], m['iph'], m['alpha'], m['d_in'])}
    rays = load_rays(out)
    rows = []
    for name, ray in sorted(rays.items()):
        if not ray['finl'] or len(ray['ems']) < 2 or name not in meta:
            continue
        iph, al, _ = meta[name]
        lphi, s, n = read_deck(os.path.join(out, 'rays', name + '.in'))
        f, sk = ray_features(lphi, s, n, al, h, r0)
        f['sin_graze'] = sin_graze(sim, itp, sim['ph'][iph], al, h, sk / r0)
        f['h'] = h
        t = ray_targets(ray['ems'], s[0], sk, r0)
        rows.append(dict(name=name, iph=iph, **f, **t))
    tab = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **tab)
    print(f"emulator table: {len(rows)} rays -> {path}")
    return tab


def neighbour_baseline(tab, train, test, key):
    """target key for the test rays from the two training columns nearest
    in phi on either side: each interpolated in alpha, then linearly in
    phi (the column index stands in for phi)"""
    cols = np.unique(tab['iph'][train])
    pred = np.full(test.sum(), np.nan)
    ti = np.nonzero(test)[0]
    for j in np.unique(tab['iph'][test]):
        lo = cols[cols < j].max() if (cols < j).any() else cols.max() - 256
        hi = cols[cols > j].min() if (cols > j).any() else cols.min() + 256
        vals = []
        for c in (lo % 256, hi % 256):
            k = train & (tab['iph'] == c) & np.isfinite(tab[key])
            o = np.argsort(tab['alpha'][k])
            vals.append((tab['alpha'][k][o], tab[key][k][o]))
        w = (j - lo) / (hi - lo)
        sel = tab['iph'][ti] == j
        a = tab['alpha'][ti[sel]]
        pred[sel] = ((1 - w) * np.interp(a, *vals[0])
                     + w * np.interp(a, *vals[1]))
    return pred


def validate(out, max_iter=400):
    """GROUP_FOLDS-fold validation by phi column of the emulator against
    the neighbour-column baseline; prints the errors and writes
    plots/emulator_validation.png"""
    from sklearn.ensemble import HistGradientBoostingRegressor
    tab = build_table(out)
    X = np.column_stack([tab[k] for k in FEATURES])
    cols = np.unique(tab['iph'])
    rng = np.random.default_rng(1)
    fold_of = dict(zip(rng.permutation(cols),
                       np.arange(len(cols)) % GROUP_FOLDS))
    fold = np.array([fold_of[c] for c in tab['iph']])
    targets = [k for k in tab if k.split('_')[0] in ('lL', 'dc', 'lw')
               and k not in FEATURES]
    pred = {k: np.full(len(X), np.nan) for k in targets}
    base = {k: np.full(len(X), np.nan) for k in targets}
    for f in range(GROUP_FOLDS):
        test, train = fold == f, fold != f
        for k in targets:
            y = tab[k]
            ok = train & np.isfinite(y)
            mdl = HistGradientBoostingRegressor(max_iter=max_iter,
                                                learning_rate=0.05,
                                                random_state=0)
            mdl.fit(X[ok], y[ok])
            pred[k][test] = mdl.predict(X[test])
            base[k][test] = neighbour_baseline(tab, train, test, k)
    print(f"{len(X)} rays, {len(cols)} phi columns, {GROUP_FOLDS} folds")
    print(f"{'target':10s} {'emulator':>28s}   {'neighbour baseline':>28s}")
    for k in targets:
        y = tab[k]
        ok = np.isfinite(y)
        e1, e2 = np.abs(pred[k] - y)[ok], np.abs(base[k] - y)[ok]
        unit = 'dex' if k.startswith('l') else 'r0'
        print(f"{k:10s} median {np.median(e1):6.3f} 84% {np.percentile(e1, 84):6.3f} {unit}"
              f"   median {np.nanmedian(e2):6.3f} 84% {np.nanpercentile(e2, 84):6.3f} {unit}")
    # what the maps and spectra sum: each held-out column's total per line
    print("per-column total luminosity (the rays' sum), |log pred/true|:")
    tot = {}
    for key in LINES:
        k = f'lL_{key}'
        y = tab[k]
        r = []
        for c in cols:
            s = (tab['iph'] == c) & np.isfinite(y)
            t = np.log10(np.sum(10 ** y[s]))
            r.append((t, np.log10(np.sum(10 ** pred[k][s])) - t,
                      np.log10(np.nansum(10 ** base[k][s])) - t))
        tot[key] = np.array(r)
        a = np.abs(tot[key][:, 1:])
        print(f"  {key:4s} emulator median {np.median(a[:, 0]):.3f} max "
              f"{a[:, 0].max():.3f} dex; baseline median "
              f"{np.median(a[:, 1]):.3f} max {a[:, 1].max():.3f} dex")
    np.savez(os.path.join(out, 'emulator', 'validation.npz'),
             **{f'pred_{k}': v for k, v in pred.items()},
             **{f'base_{k}': v for k, v in base.items()})
    plt = plot_style()
    fig, axs = plt.subplots(2, 3, figsize=(20, 13))
    for ax, key in zip(axs.ravel(), LINES):
        k = f'lL_{key}'
        y = tab[k]
        ax.plot(y, base[k], '.', color='0.6', ms=3, label=r'\textbf{neighbours}')
        ax.plot(y, pred[k], '.', color='orangered', ms=3,
                label=r'\textbf{emulator}')
        lo, hi = np.nanpercentile(y, [0.5, 99.9])
        ax.plot([lo, hi], [lo, hi], 'k-', lw=1)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_title(r'\textbf{' + LINE_TEX[key] + r'}', fontsize=23)
        ax.set_xlabel(r'$\log L_{\rm Cloudy}$')
        ax.tick_params(which='both', direction='in', labelsize=17)
    for ax in axs[:, 0]:
        ax.set_ylabel(r'$\log L_{\rm predicted}$')
    axs[0, 0].legend(fontsize=17, frameon=False, loc='upper left',
                     markerscale=4)
    ax = axs[1, 2]
    for key, c in zip(LINES, ('k', 'gray', 'dodgerblue', 'orangered',
                              'darkviolet')):
        ax.plot(tot[key][:, 0] - tot[key][:, 0].max(), tot[key][:, 1], 'o',
                color=c, ms=5, label=r'\textbf{' + LINE_TEX[key] + r'}')
    ax.axhline(0, color='k', lw=1)
    ax.set_xlabel(r'$\log L_{\rm column} - \max$')
    ax.set_ylabel(r'\textbf{column total: emulator} $-$ \textbf{Cloudy [dex]}')
    ax.legend(fontsize=15, frameon=False, ncol=2)
    ax.tick_params(which='both', direction='in', labelsize=17)
    fig.tight_layout()
    pdir = os.path.join(out, 'plots')
    os.makedirs(pdir, exist_ok=True)
    path = os.path.join(pdir, 'emulator_validation.png')
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {path}")


def validate_heights(runs=HEIGHT_RUNS, max_iter=400):
    """Leave one lamp height out: for each run in runs (one lamp height
    each), train on all the others' rays and predict its rays; prints the
    errors per target and of the run's total luminosity per line (what its
    maps and spectra sum, weighted by each ray's solid angle)"""
    from sklearn.ensemble import HistGradientBoostingRegressor
    runs = [r for r in runs if os.path.exists(os.path.join(r, 'rays_meta.npz'))]
    tabs = [build_table(r) for r in runs]
    doms = []
    for r, t in zip(runs, tabs):
        m = np.load(os.path.join(r, 'rays_meta.npz'))
        dom = dict(zip(m['name'].astype(str), m['domega']))
        doms.append(np.array([dom[str(n)] for n in t['name']]))
    targets = [k for k in tabs[0] if k.split('_')[0] in ('lL', 'dc', 'lw')
               and k not in FEATURES]
    print(f"leave-one-height-out over {len(runs)} runs: " + ', '.join(
        f"h = {float(t['h'][0]):.2f} r0 ({len(t['h'])} rays)" for t in tabs))
    for i, (r, t) in enumerate(zip(runs, tabs)):
        rest = [tt for j, tt in enumerate(tabs) if j != i]
        X = np.column_stack([np.concatenate([tt[k] for tt in rest])
                             for k in FEATURES])
        Xt = np.column_stack([t[k] for k in FEATURES])
        lo, hi = np.percentile(X, 1, axis=0), np.percentile(X, 99, axis=0)
        cov = [j for j, k in enumerate(FEATURES) if k != 'h']  # h always is
        outside = np.mean(((Xt < lo) | (Xt > hi))[:, cov].any(axis=1))
        print(f"== held out h = {float(t['h'][0]):.2f} r0 ({r}); "
              f"rays with a feature outside the others' 1-99%: "
              f"{outside:.0%}")
        for k in targets:
            y = np.concatenate([tt[k] for tt in rest])
            ok = np.isfinite(y)
            mdl = HistGradientBoostingRegressor(max_iter=max_iter,
                                                learning_rate=0.05,
                                                random_state=0)
            mdl.fit(X[ok], y[ok])
            p = mdl.predict(Xt)
            e = np.abs(p - t[k])[np.isfinite(t[k])]
            line = (f"  {k:8s} median {np.median(e):6.3f} 84% "
                    f"{np.percentile(e, 84):6.3f} "
                    f"{'dex' if k.startswith('l') else 'r0'}")
            if k.startswith('lL_'):
                s_ = np.isfinite(t[k])
                tot = np.log10(np.sum(10 ** t[k][s_] * doms[i][s_]))
                pt = np.log10(np.sum(10 ** p[s_] * doms[i][s_]))
                line += f";  run total {pt - tot:+.3f} dex"
            print(line, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('step', choices=['table', 'validate', 'heights'])
    ap.add_argument('--out', default=None,
                    help='table / validate: the run directory')
    ap.add_argument('--runs', nargs='+', default=list(HEIGHT_RUNS),
                    help='heights: one run per lamp height')
    ap.add_argument('--max-iter', type=int, default=400,
                    help='boosting iterations per target model')
    a = ap.parse_args()
    if a.step == 'table':
        build_table(a.out)
    elif a.step == 'validate':
        validate(a.out, a.max_iter)
    else:
        validate_heights(a.runs, a.max_iter)


if __name__ == '__main__':
    main()
