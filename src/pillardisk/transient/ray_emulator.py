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
  lands      1 if n_H ever exceeds N_SURF (else the skin is the n_H peak)
  lNS<f>     log Stromgren column: the column from the entry within which
             case-B recombinations (alpha_B 1.1 n_H^2 s^2 ds) use up the
             fraction f/100 of the entry's ionising photons (phi s0^2;
             SFRAC: 100% to 1%); log column at the end + 0.5 if never
             reached (whether and where the front breaks through the rim)
  lrec       log of the ray's total recombinations / ionising photons
Only what Cloudy's 1D slab sees: no lamp height, landing position or
incidence angle (the geometry enters afterwards, placing and projecting the
layer), so a new lamp height is emulated wherever its rays fall inside the
training rays' feature space.
Targets (from the run's ray cache), per line (Ha, Hb, C IV, Mg II, He II),
slab quantities:
  lI_<line>  log sum 4 pi j (s/s_skin)^2 ds [erg/s/cm^2]: per unit area at
             the skin (Cloudy dilutes spherically from the entry); the ray's
             log sum 4 pi j s^2 ds = lI_<line> + 2 log s_skin exactly
  lq<f>_<line>  log of the column from the entry [cm^-2] within which the
             fraction f/1000 of I is emitted (QFRAC: 0.5% to 99.5%): the
             emission's whole distribution in column, tails included (at a
             low lamp a ray lights the rim, crosses thin gas and lights a
             second layer far out: one emission curve in column, two
             layers in distance)
  lTe_<line> log T_e at the emission's (depth) centroid [K]
  lrat_<line>  log ratio of the doublet's components (C IV, Mg II)
and for C IV's escape (as cloudy_rays.civ_layers, from n(C3+)):
  lnl_c4     log n(C3+) at the centroid [cm^-3]
  lNab_c4    log C3+ column from the entry to the centroid [cm^-2]
  lNtot_c4   log the ray's whole C3+ column
  lthk_c4    log the path where the C3+ fraction exceeds 0.5 [cm]
(positions in s follow from a ray's own column profile, its deck)

`emulate` trains on runs and writes an emulated run: every ray of another
run's plan (its decks; Cloudy need not have run them) gets a synthetic
emissivity profile per line (a two-sided Gaussian in column through the
predicted column quantiles: a monotone spline of the emitted fraction
against log column, mapped to depth through the deck's own density
profile, normalised to I) and its emulated
layer; cloudy_rays' map, escape and profile steps then run on it
(<out>/emulated.npz, read by load_rays and the layer functions).

`validate` holds out whole phi columns (GROUP_FOLDS folds; neighbouring rays
in a column are correlated) and compares, ray by ray and per held-out column
total, the emulator (one gradient-boosted tree model per target) with a
baseline that interpolates the two neighbouring training columns in alpha
and phi. `heights` trains on all but one lamp height (HEIGHT_RUNS) and
predicts the one left out: can a height be emulated from the others?

Run:  python3 -m pillardisk.transient.ray_emulator validate --out RUN_DIR
      python3 -m pillardisk.transient.ray_emulator heights
      python3 -m pillardisk.transient.ray_emulator heightsplot   (replot, no training)
      python3 -m pillardisk.transient.ray_emulator emulate --src RUN --out DIR
      python3 -m pillardisk.transient.ray_emulator compare --out DIR   (face-on maps)
      python3 -m pillardisk.transient.ray_emulator dense --h 0.5 --nrays 30000  (+ map)
      python3 -m pillardisk.transient.ray_emulator redo --run h0.5_n30000   (re-emulate)
      python3 -m pillardisk.transient.ray_emulator synth --run NAME  (rebuild + map, no training)
      python3 -m pillardisk.transient.ray_emulator gapfill --src DENSE_RUN   (Cloudy where
          the emulator is blind; picked up by the next emulate)
(the tree models use EMU_THREADS threads: all cores fight on a busy machine)
"""

import argparse
import glob
import os
import shutil

EMU_THREADS = 8
os.environ.setdefault('OMP_NUM_THREADS', str(EMU_THREADS))

import numpy as np                                        # noqa: E402

from pillardisk.transient.cloudy_rays import (load_rays, N_SURF, plot_style,  # noqa
                                   CACHE_OVR, NEMS)
from pillardisk.transient.sim_spectrum import CONFIG, LD_CM                   # noqa

LINES = dict(ha=(1,), hb=(2,), c4=(3, 4), mg2=(5, 6), he2=(7,))
LINE_TEX = dict(ha=r'H$\bm{\alpha}$', hb=r'H$\bm{\beta}$', c4=r'C\,IV',
                mg2=r'Mg\,II', he2=r'He\,II\,1640')
FEATURES = ('lphi_skin', 'lN_pre', 'lL_skin', 'ln_p1', 'ln_p3', 'lN_post',
            'lnmax', 'lands')
SFRAC = (1.0, 0.3, 0.1, 0.03, 0.01)    # photon fractions at the lNS features
ALPHA_B = 2.6e-13                       # case-B H recombination [cm^3/s]
FEATURES += tuple(f'lNS{round(f * 100):03d}' for f in SFRAC) + ('lrec',)
QFRAC = (0.005, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.98, 0.99,
         0.995)             # emitted fractions at the lq targets
Q_T = tuple(f'lq{round(q * 1000):03d}' for q in QFRAC)
LINE_T = ('lI',) + Q_T + ('lTe',)                  # every line
C4_T = ('lnl', 'lNab', 'lNtot', 'lthk')            # C IV's C3+ layer
TARGETS = LINE_T + ('lrat',) + C4_T
TABLE_VERSION = 5           # + Stromgren-column features (SFRAC)
RUNS = '/data2/jhuang/runs/cloudy_rays'
HEIGHT_RUNS = (f'{RUNS}/train/h0.0', f'{RUNS}/train/h0.2',
               f'{RUNS}/edd10_h0.5_adapt', f'{RUNS}/train/h1.0')
VAL_RUNS = {'h0.0': HEIGHT_RUNS[0], 'h0.2': HEIGHT_RUNS[1],
            'h0.5': HEIGHT_RUNS[2], 'h1.0': HEIGHT_RUNS[3]}
GROUP_FOLDS = 8
MIN_RAYS = 50               # heights: a run needs this many finished rays
HEIGHTS_NPZ = f'{RUNS}/train/emulator_heights.npz'   # heights' predictions


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


def column(s, n):
    """column [cm^-2] from the entry along the deck's (s, n_H)"""
    return np.concatenate([[0], np.cumsum(0.5 * (n[1:] + n[:-1])
                                          * np.diff(s))])


def ray_features(lphi, s, n, r0):
    """the FEATURES of one ray (see the module docstring), and s_skin"""
    col = column(s, n)
    hit = np.nonzero(n >= N_SURF)[0]
    k = hit[0] if len(hit) else int(np.argmax(n))
    sk = s[k]
    lnn = np.log(n)
    j0, j1 = max(k - 3, 0), min(k + 3, len(s) - 1)
    g = (lnn[j1] - lnn[j0]) / max(s[j1] - s[j0], 1.0)  # d ln n / ds [1/cm]
    lsc = np.log10(np.clip(1 / max(abs(g), 1e-30) / r0, 1e-6, 10))
    at = lambda ds: np.log10(np.interp(sk + ds * r0, s, n))
    g2 = ALPHA_B * 1.1 * n ** 2 * s ** 2
    rec = np.concatenate([[0], np.cumsum(0.5 * (g2[1:] + g2[:-1])
                                         * np.diff(s))])
    need = 10 ** lphi * s[0] ** 2
    strom = {f'lNS{round(f * 100):03d}':
             (np.log10(np.interp(f * need, rec, col) + 1e10)
              if rec[-1] > f * need else np.log10(col[-1] + 1e10) + 0.5)
             for f in SFRAC}
    return dict(**strom, lrec=np.log10(max(rec[-1], 1e-30) / need),
                lphi_skin=lphi + 2 * np.log10(s[0] / sk),
                lN_pre=np.log10(col[k] + 1e10),
                lL_skin=lsc, ln_p1=at(0.01), ln_p3=at(0.03),
                lN_post=np.log10(col[-1] - col[k] + 1e10),
                lnmax=np.log10(n.max()),
                lands=float(len(hit) > 0)), sk


def ray_targets(ray, s_deck, n_deck, sk):
    """the TARGETS of one ray from its cache entry (.ems: depth, 4 pi j
    columns; .ovr; C/H) and its deck's (s, n_H)"""
    e, ov = ray['ems'], ray['ovr']
    s = s_deck[0] + e[:, 0]
    ds = np.gradient(e[:, 0])
    col = np.interp(s, s_deck, column(s_deck, n_deck))
    oc = {k: ov[:, 1 + i] for i, k in enumerate(CACHE_OVR)}
    out = {}
    for key, cols in LINES.items():
        w = e[:, list(cols)].sum(axis=1) * (s / sk) ** 2 * ds
        tot = w.sum()
        keys = LINE_T + (('lrat',) if len(cols) == 2 else ()) + \
            (C4_T if key == 'c4' else ())
        out.update({f'{t}_{key}': np.nan for t in keys})
        if not tot > 0:
            continue
        cw = np.cumsum(w) / tot
        dc = np.sum(w * e[:, 0]) / tot
        out[f'lI_{key}'] = np.log10(tot)
        for t, q in zip(Q_T, QFRAC):
            out[f'{t}_{key}'] = np.log10(max(np.interp(q, cw, col), 1e10))
        out[f'lTe_{key}'] = np.log10(np.interp(dc, ov[:, 0], oc['Te']))
        if len(cols) == 2:
            a, b = (np.sum(e[:, c] * ds) for c in cols)
            if a > 0 and b > 0:
                out[f'lrat_{key}'] = np.log10(a / b)
        if key == 'c4':
            nl = oc['C4'] * oc['hden'] * ray['cha']
            dd = np.gradient(ov[:, 0])
            Ncum = np.cumsum(nl * dd)
            out['lnl_c4'] = np.log10(max(np.interp(dc, ov[:, 0], nl),
                                         1e-10))
            out['lNab_c4'] = np.log10(max(np.interp(dc, ov[:, 0], Ncum),
                                          1.0))
            out['lNtot_c4'] = np.log10(max(Ncum[-1], 1.0))
            out['lthk_c4'] = np.log10(max(np.sum(dd[oc['C4'] > 0.5]), 1.0))
    return out


def build_table(out):
    """features and targets of every finished ray of the run (adaptive
    layout: rays_meta.npz), cached as <out>/emulator/table.npz"""
    path = os.path.join(out, 'emulator', 'table.npz')
    rays = load_rays(out)
    nfin = sum(bool(r['finl']) for r in rays.values())
    if os.path.exists(path):
        d = np.load(path)
        # rebuilt when more rays have finished (a run still going)
        if 'version' in d.files and int(d['version']) == TABLE_VERSION \
                and int(d['nfin']) == nfin:
            return {k: d[k] for k in d.files
                    if k not in ('nfin', 'version')}
    if nfin == 0:
        return None
    m = np.load(os.path.join(out, 'rays_meta.npz'))
    h, r0 = float(m['hlamp_r0']), CONFIG['R0_LD'] * LD_CM
    meta = {str(n): (int(p), float(a), float(d0)) for n, p, a, d0 in
            zip(m['name'], m['iph'], m['alpha'], m['d_in'])}
    rows = []
    for name, ray in sorted(rays.items()):
        if not ray['finl'] or len(ray['ems']) < 2 or name not in meta:
            continue
        iph, al, _ = meta[name]
        lphi, s, n = read_deck(os.path.join(out, 'rays', name + '.in'))
        f, sk = ray_features(lphi, s, n, r0)
        t = ray_targets(ray, s, n, sk)
        # alpha (the baseline), h and s_skin (L = s_skin^2 I) are labels,
        # not features
        rows.append(dict(name=name, iph=iph, alpha=al, h=h,
                         lsk=np.log10(sk), **f, **t))
    tab = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, nfin=nfin, version=TABLE_VERSION, **tab)
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
    targets = [k for k in tab if k.split('_')[0] in TARGETS]
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
        print(f"{k:10s} median {np.median(e1):6.3f} 84% {np.percentile(e1, 84):6.3f} dex"
              f"   median {np.nanmedian(e2):6.3f} 84% {np.nanpercentile(e2, 84):6.3f} dex")
    # what the maps and spectra sum: each held-out column's total per line
    print("per-column total luminosity (the rays' sum), |log pred/true|:")
    tot = {}
    for key in LINES:
        k = f'lI_{key}'
        y = tab[k] + 2 * tab['lsk']                      # log L = lI + 2 lsk
        yp, yb = pred[k] + 2 * tab['lsk'], base[k] + 2 * tab['lsk']
        r = []
        for c in cols:
            s = (tab['iph'] == c) & np.isfinite(y)
            t = np.log10(np.sum(10 ** y[s]))
            r.append((t, np.log10(np.sum(10 ** yp[s])) - t,
                      np.log10(np.nansum(10 ** yb[s])) - t))
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
        k = f'lI_{key}'
        y = tab[k]
        ax.plot(y, base[k], '.', color='0.6', ms=3, label=r'\textbf{neighbours}')
        ax.plot(y, pred[k], '.', color='orangered', ms=3,
                label=r'\textbf{emulator}')
        lo, hi = np.nanpercentile(y, [0.5, 99.9])
        ax.plot([lo, hi], [lo, hi], 'k-', lw=1)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_title(r'\textbf{' + LINE_TEX[key] + r'}', fontsize=23)
        ax.set_xlabel(r'$\log I_{\rm Cloudy}$')
        ax.tick_params(which='both', direction='in', labelsize=17)
    for ax in axs[:, 0]:
        ax.set_ylabel(r'$\log I_{\rm predicted}$')
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
    maps and spectra sum, weighted by each ray's solid angle); saves every
    held-out ray's prediction to HEIGHTS_NPZ and plots them (heights_plot)"""
    from sklearn.ensemble import HistGradientBoostingRegressor
    runs = [r for r in runs if os.path.exists(os.path.join(r, 'rays_meta.npz'))]
    tabs = [build_table(r) for r in runs]
    # heights still running count with the rays finished so far
    keep = [t is not None and len(t['h']) >= MIN_RAYS for t in tabs]
    for r, k in zip(runs, keep):
        if not k:
            print(f"skipped {r}: fewer than {MIN_RAYS} finished rays")
    runs = [r for r, k in zip(runs, keep) if k]
    tabs = [t for t, k in zip(tabs, keep) if k]
    doms, save = [], {}
    for r, t in zip(runs, tabs):
        m = np.load(os.path.join(r, 'rays_meta.npz'))
        dom = dict(zip(m['name'].astype(str), m['domega']))
        doms.append(np.array([dom[str(n)] for n in t['name']]))
    targets = [k for k in tabs[0] if k.split('_')[0] in TARGETS]
    print(f"leave-one-height-out over {len(runs)} runs: " + ', '.join(
        f"h = {float(t['h'][0]):.2f} r0 ({len(t['h'])} rays)" for t in tabs))
    for i, (r, t) in enumerate(zip(runs, tabs)):
        rest = [tt for j, tt in enumerate(tabs) if j != i]
        X = np.column_stack([np.concatenate([tt[k] for tt in rest])
                             for k in FEATURES])
        Xt = np.column_stack([t[k] for k in FEATURES])
        lo, hi = np.percentile(X, 1, axis=0), np.percentile(X, 99, axis=0)
        out_ = ((Xt < lo) | (Xt > hi)).any(axis=1)
        outside = np.mean(out_)
        save[f'r{i}_h'] = float(t['h'][0])
        save[f'r{i}_outside'] = out_
        save[f'r{i}_lsk'] = t['lsk']
        save[f'r{i}_domega'] = doms[i]
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
            save[f'r{i}_true_{k}'], save[f'r{i}_pred_{k}'] = t[k], p
            e = np.abs(p - t[k])[np.isfinite(t[k])]
            line = (f"  {k:8s} median {np.median(e):6.3f} 84% "
                    f"{np.percentile(e, 84):6.3f} dex")
            if k.startswith('lI_'):
                s_ = np.isfinite(t[k])
                w = doms[i][s_] * 10 ** (2 * t['lsk'][s_])  # L = s_skin^2 I
                tot = np.log10(np.sum(10 ** t[k][s_] * w))
                pt = np.log10(np.sum(10 ** p[s_] * w))
                line += f";  run total {pt - tot:+.3f} dex"
            print(line, flush=True)
    np.savez(HEIGHTS_NPZ, nrun=len(runs), runs=np.array(runs), **save)
    print(f"wrote {HEIGHTS_NPZ}")
    heights_plot()


def heights_plot(path=HEIGHTS_NPZ):
    """the held-out rays' emulated against Cloudy I, lines (columns) x
    held-out heights (rows); rays with a feature outside the training
    runs' 1-99% in grey; each panel's median per-ray error and run-total
    error. Writes plots/emulator_heights.png next to HEIGHTS_NPZ"""
    d = np.load(path)
    nrun = int(d['nrun'])
    order = sorted(range(nrun), key=lambda i: float(d[f'r{i}_h']))
    plt = plot_style()
    fig, axs = plt.subplots(nrun, len(LINES), figsize=(5 * len(LINES),
                                                       5 * nrun),
                            sharex='col', sharey='col')
    axs = np.atleast_2d(axs)
    for c, key in enumerate(LINES):
        k = f'lI_{key}'
        allv = np.concatenate([d[f'r{i}_true_{k}'] for i in order])
        lo, hi = np.nanpercentile(allv, [0.5, 99.9])
        lo, hi = lo - 0.2, hi + 0.2
        for row, i in enumerate(order):
            ax = axs[row, c]
            y, p = d[f'r{i}_true_{k}'], d[f'r{i}_pred_{k}']
            ok = np.isfinite(y)
            o = d[f'r{i}_outside'].astype(bool)
            ax.plot(y[ok & ~o], p[ok & ~o], '.', color='dodgerblue', ms=4,
                    label=r'\textbf{inside}')
            ax.plot(y[ok & o], p[ok & o], '.', color='0.6', ms=4,
                    label=r'\textbf{outside}')
            ax.plot([lo, hi], [lo, hi], 'k-', lw=1)
            w = d[f'r{i}_domega'][ok] * 10 ** (2 * d[f'r{i}_lsk'][ok])
            tot = (np.log10(np.sum(10 ** p[ok] * w))
                   - np.log10(np.sum(10 ** y[ok] * w)))
            ax.text(0.04, 0.96, r'\textbf{median %.3f dex}' % np.median(
                np.abs(p - y)[ok]) + '\n' + r'\textbf{total %+.3f dex}'
                    % tot, transform=ax.transAxes, va='top', fontsize=17)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.minorticks_on()
            ax.tick_params(which='both', direction='in', labelsize=17,
                           top=True, right=True)
            ax.tick_params(which='major', length=10, width=2)
            ax.tick_params(which='minor', length=5, width=1.5)
            if row == 0:
                ax.set_title(r'\textbf{' + LINE_TEX[key] + r'}', fontsize=23)
            if c == len(LINES) - 1:
                ax.text(1.04, 0.5, r'$\bm{h = %.1f}$\,\textbf{ld}'
                        % (float(d[f'r{i}_h']) * CONFIG['R0_LD']),
                        transform=ax.transAxes, rotation=270, va='center',
                        fontsize=21)
        axs[-1, c].set_xlabel(r'$\log I_{\rm Cloudy}$')
    for ax in axs[:, 0]:
        ax.set_ylabel(r'$\log I_{\rm emulated}$')
    axs[0, 0].legend(fontsize=15, frameon=False, loc='lower right',
                     markerscale=3)
    fig.tight_layout(h_pad=0.3, w_pad=0.3)
    pdir = os.path.join(os.path.dirname(path), 'plots')
    os.makedirs(pdir, exist_ok=True)
    out = os.path.join(pdir, 'emulator_heights.png')
    fig.savefig(out, dpi=100, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {out}")


def fit_models(tabs, max_iter=400):
    """{target: model} trained on the rays of tables tabs"""
    from sklearn.ensemble import HistGradientBoostingRegressor
    X = np.column_stack([np.concatenate([t[k] for t in tabs])
                         for k in FEATURES])
    models = {}
    for k in [k for k in tabs[0] if k.split('_')[0] in TARGETS]:
        y = np.concatenate([t[k] for t in tabs])
        ok = np.isfinite(y)
        models[k] = HistGradientBoostingRegressor(
            max_iter=max_iter, learning_rate=0.05, random_state=0).fit(
                X[ok], y[ok])
    return models


SYN_DS = 0.005              # [r0] largest step of a synthetic profile (the
                            # sim's radial cells are 0.012 r0)


def column_cdf(lq):
    """the emitted fraction F(log N) from one line's log column quantiles
    lq (at QFRAC; made non-decreasing): a monotone (PCHIP) spline through
    them, F = 0 and 1 one quantile spacing beyond the ends. Returns
    (F, its log N range)"""
    from scipy.interpolate import PchipInterpolator
    x = np.maximum.accumulate(np.asarray(lq, float))
    x = x + 1e-4 * np.arange(len(x))          # strictly increasing
    x = np.concatenate([[x[0] - max(x[1] - x[0], 0.01)], x,
                        [x[-1] + max(x[-1] - x[-2], 0.01)]])
    f = np.concatenate([[0.0], QFRAC, [1.0]])
    return PchipInterpolator(x, f, extrapolate=False), (x[0], x[-1])


def synth_ray(p, s, n, sk, npg=60):
    """one emulated ray from its predicted targets p {target: value} and
    its deck (s, n_H; s_skin): (ems (npt, NEMS) on a depth grid, {line:
    (d_centroid, d_5%, d_95%) [cm]}). Each line: the emitted fraction
    against column (column_cdf of its lq quantiles), differentiated and
    placed along the ray by the deck's column profile, normalised to sum
    4 pi j (s/s_skin)^2 ds = I. Low density stretches (little column) get
    little emission, so a curve continuous in column can split into
    separate layers in depth."""
    col = column(s, n)
    nmax = col[-1]
    cdf, grids = {}, [[0.0, nmax]]
    ds = SYN_DS * CONFIG['R0_LD'] * LD_CM
    for key in LINES:
        F, (xa, xb) = column_cdf([p[f'{t}_{key}'] for t in Q_T])
        cdf[key] = F
        na, nb = np.clip([10 ** xa, 10 ** xb], 0, nmax)
        # log-even in column, and at most SYN_DS apart in s: the face-on
        # map puts each point in one radial cell, so a wide (grazing)
        # layer sampled coarser than the cells leaves empty cells
        grids.append(np.logspace(np.log10(max(na, 1e10)),
                                 np.log10(max(nb, 1e10)), npg))
        sa, sb = np.interp([na, nb], col, s)
        grids.append(np.interp(np.arange(sa, sb, ds), s, col))
        grids.append(col[(col >= na) & (col <= nb)])
    Ng = np.unique(np.clip(np.concatenate(grids), 0, nmax))
    sg = np.interp(Ng, col, s)
    d = sg - s[0]
    dd = np.gradient(d)
    ng = np.interp(sg, s, n)
    lN = np.log10(np.maximum(Ng, 1e-30))
    e = np.zeros((len(d), NEMS))
    e[:, 0] = d
    lay = {}
    for key, cols in LINES.items():
        F = cdf[key]
        # dF/dN = dF/dlogN / (N ln 10); zero outside the spline's range
        q = np.nan_to_num(F.derivative()(lN)) / (np.maximum(Ng, 1e-30)
                                                 * np.log(10))
        q = np.maximum(q, 0)
        w = q * ng * dd                       # = j (s/s_skin)^2 ds
        tot = w.sum()
        j = q * ng * (sk / sg) ** 2 * 10 ** p[f'lI_{key}'] / max(tot,
                                                                 1e-300)
        if len(cols) == 2:
            r = 10 ** p[f'lrat_{key}']
            e[:, cols[0]], e[:, cols[1]] = j * r / (1 + r), j / (1 + r)
        else:
            e[:, cols[0]] = j
        dc = np.sum(w * d) / max(tot, 1e-300)
        cw = np.cumsum(w) / max(tot, 1e-300)
        lay[key] = (dc,) + tuple(np.interp([0.05, 0.95], cw, d))
    return e, lay


SMOOTH_N = 5                # rays in the running median along each column


MODELS_DIR = f'{RUNS}/train/emulator_models'     # trained models, reused
GAP_N = 300                 # gapfill: Cloudy rays per gap set
GAP_NCPU = 40               # gapfill: cores (niced)


def training_runs(exclude=None):
    """the emulator's training runs: HEIGHT_RUNS and every gapfill set
    (<RUNS>/train/gap_*), less the run exclude"""
    gaps = sorted(os.path.dirname(p) for p in
                  glob.glob(f'{RUNS}/train/gap_*/rays_meta.npz'))
    ex = os.path.abspath(exclude) if exclude else None
    return [r for r in list(HEIGHT_RUNS) + gaps if os.path.abspath(r) != ex]


def gapfill(src, n=GAP_N, out=None, ncpu=GAP_NCPU, run=True):
    """A training set where the emulator is blind: of the emulated run
    src's rays whose features fall outside the training rays' 1-99%
    (emulated_pred.npz 'outside'), n spread over that region (farthest-
    point sampling in features standardised by the training rays), their
    decks copied to out (default <RUNS>/train/gap_<src name>) with their
    rays_meta rows; then Cloudy on them (niced, ncpu at a time,
    unfinished only: a stopped run resumes). training_runs() picks the set
    up, so the next emulate retrains with it."""
    from pillardisk.transient.cloudy_rays import run_rays
    src = os.path.abspath(src)
    out = out or f'{RUNS}/train/gap_{os.path.basename(src)}'
    if not os.path.exists(os.path.join(out, 'rays_meta.npz')):
        d = np.load(os.path.join(src, 'emulated_pred.npz'))
        names, X = d['name'].astype(str), d['X']
        cand = np.nonzero(d['outside'])[0]
        Xt = np.column_stack([np.concatenate(
            [build_table(r)[k] for r in training_runs(exclude=src)])
            for k in FEATURES])
        mu, sd = Xt.mean(axis=0), np.maximum(Xt.std(axis=0), 1e-6)
        Z = (X[cand] - mu) / sd
        pick = [int(np.argmax(np.linalg.norm(Z, axis=1)))]
        dist = np.linalg.norm(Z - Z[pick[0]], axis=1)
        while len(pick) < min(n, len(cand)):
            pick.append(int(np.argmax(dist)))
            dist = np.minimum(dist, np.linalg.norm(Z - Z[pick[-1]], axis=1))
        chosen = set(names[cand[pick]])
        m = np.load(os.path.join(src, 'rays_meta.npz'))
        keep = np.array([str(nm) in chosen for nm in m['name']])
        os.makedirs(os.path.join(out, 'rays'), exist_ok=True)
        for nm in m['name'][keep]:
            shutil.copy(os.path.join(src, 'rays', f'{nm}.in'),
                        os.path.join(out, 'rays', f'{nm}.in'))
        np.savez(os.path.join(out, 'rays_meta.npz'),
                 **{k: (m[k][keep] if m[k].ndim else m[k]) for k in m.files})
        print(f"wrote {keep.sum()} decks to {out}: of {len(cand)} rays of "
              f"{src} outside the training range ({d['outside'].mean():.0%})",
              flush=True)
    if run:
        os.nice(19)
        run_rays(out, ncpu)


def trained_models(train, max_iter=400):
    """{target: model} trained on the runs train, loaded from MODELS_DIR
    when a saved set matches (same runs, rays, table version and max_iter),
    else fitted (fit_models) and saved there (one file per training set)"""
    import hashlib
    import joblib
    tabs = [build_table(r) for r in train]
    key = repr((sorted(os.path.abspath(r) for r in train),
                [len(t['name']) for t in tabs], TABLE_VERSION, max_iter))
    path = os.path.join(MODELS_DIR, hashlib.md5(key.encode()).hexdigest()[:12]
                        + '.joblib')
    if os.path.exists(path):
        print(f"loaded models {path}", flush=True)
        return joblib.load(path)
    print("training on " + ', '.join(
        f"h = {float(t['h'][0]):.2f} r0 ({len(t['h'])} rays)" for t in tabs),
        flush=True)
    models = fit_models(tabs, max_iter)
    os.makedirs(MODELS_DIR, exist_ok=True)
    joblib.dump(models, path)
    print(f"saved models {path}", flush=True)
    return models


def emulate(src, out, train=None, max_iter=400):
    """An emulated run in out: the ray plan of run src (rays_meta.npz and
    its decks), every ray's targets predicted by models trained on the
    runs train (default training_runs() without src; trained_models: fitted
    once per training set, then reused), saved raw to
    <out>/emulated_pred.npz (with each ray's s_skin); then synth"""
    src = os.path.abspath(src)
    if train is None:
        train = training_runs(exclude=src)
    models = trained_models(train, max_iter)
    m = np.load(os.path.join(src, 'rays_meta.npz'))
    r0 = CONFIG['R0_LD'] * LD_CM
    names = sorted(m['name'].astype(str))
    feats, sks = [], []
    for name in names:
        lphi, s, n = read_deck(os.path.join(src, 'rays', name + '.in'))
        f, sk = ray_features(lphi, s, n, r0)
        feats.append([f[k] for k in FEATURES])
        sks.append(sk)
    X = np.array(feats)
    Xt = np.column_stack([np.concatenate([build_table(r)[k] for r in train])
                          for k in FEATURES])
    lo, hi = np.percentile(Xt, 1, axis=0), np.percentile(Xt, 99, axis=0)
    outside = ((X < lo) | (X > hi)).any(axis=1)
    print(f"rays with a feature outside the training rays' 1-99%: "
          f"{outside.mean():.0%} (trust drops there: trees do not "
          f"extrapolate)", flush=True)
    cha = next(r['cha'] for r in load_rays(train[0]).values()
               if np.isfinite(r['cha']))
    os.makedirs(out, exist_ok=True)
    if os.path.abspath(out) != src:
        shutil.copy(os.path.join(src, 'rays_meta.npz'),
                    os.path.join(out, 'rays_meta.npz'))
    path = os.path.join(out, 'emulated_pred.npz')
    np.savez(path, name=np.array(names), src=src, train=np.array(train),
             version=TABLE_VERSION,
             cha=cha, sk=np.array(sks), X=X, outside=outside,
             **{f'pred_{k}': mdl.predict(X) for k, mdl in models.items()})
    print(f"wrote {path}: {len(names)} rays' predicted targets", flush=True)
    synth(out)


def smooth_columns(names, meta, v, nwin=SMOOTH_N):
    """v (per ray, aligned with names) with a running median over nwin
    neighbouring rays in alpha within each phi column (ends: the nearest
    rays); the trees' piecewise-constant predictions step between
    neighbours, Cloudy's own outputs do not"""
    from scipy.ndimage import median_filter
    al = dict(zip(meta['name'].astype(str), meta['alpha']))
    iph = np.array([int(n[5:8]) for n in names])
    a = np.array([al[n] for n in names])
    out = v.copy()
    for c in np.unique(iph):
        k = np.nonzero(iph == c)[0]
        k = k[np.argsort(a[k])]
        out[k] = median_filter(v[k], size=min(nwin, len(k)), mode='nearest')
    return out


def synth(out, nwin=SMOOTH_N):
    """<out>/emulated.npz from the raw predictions <out>/emulated_pred.npz
    (no training): each target smoothed along its phi column
    (smooth_columns, nwin rays; 1 = none), then each ray's synthetic .ems
    (synth_ray; depths float64, emissivities float32) and, per line, its
    layer (d_c, d_lo, d_hi [cm], T_e) and C IV's C3+ quantities, which
    cloudy_rays reads in place of Cloudy's output"""
    d = np.load(os.path.join(out, 'emulated_pred.npz'))
    src = str(d['src'])
    names = d['name'].astype(str)
    meta = np.load(os.path.join(out, 'rays_meta.npz'))
    pred = {k[5:]: d[k] for k in d.files if k.startswith('pred_')}
    if nwin > 1:
        pred = {k: smooth_columns(names, meta, v, nwin)
                for k, v in pred.items()}
    ems, lay = [], {key: [] for key in LINES}
    for i, name in enumerate(names):
        _, s, n = read_deck(os.path.join(src, 'rays', name + '.in'))
        e, ly = synth_ray({k: v[i] for k, v in pred.items()}, s, n,
                          d['sk'][i])
        ems.append(e)
        for key in LINES:
            lay[key].append(ly[key])
    E = np.concatenate(ems)
    save = dict(name=names, src=src, train=d['train'], cha=d['cha'],
                nwin=nwin, n_ems=np.array([len(e) for e in ems]),
                ems_d=E[:, 0], ems_v=E[:, 1:].astype(np.float32))
    for key in LINES:
        a = np.array(lay[key])
        save.update({f'dc_{key}': a[:, 0], f'dlo_{key}': a[:, 1],
                     f'dhi_{key}': a[:, 2],
                     f'Te_{key}': 10 ** pred[f'lTe_{key}']})
    save.update(nl_c4=10 ** pred['lnl_c4'], Nabove_c4=10 ** pred['lNab_c4'],
                Ntot_c4=10 ** pred['lNtot_c4'],
                thick_c4=10 ** pred['lthk_c4'])
    path = os.path.join(out, 'emulated.npz')
    np.savez(path, **save)
    print(f"wrote {path}: {len(names)} emulated rays of {src} (running "
          f"median over {nwin} rays; {os.path.getsize(path) / 1e6:.0f} MB)")


def dense(h_ld, nrays, nphi, out=None, max_iter=400):
    """A dense emulated run: cloudy_rays' adaptive plan of ~nrays rays in
    nphi phi columns (at most one per sim phi cell) for a lamp at h_ld
    [ld], traced (decks only, no Cloudy) into out, then emulated there by
    models trained on all training_runs(). Returns out."""
    from pillardisk.transient.cloudy_rays import write_adaptive
    out = out or f'{RUNS}/emulated/h{h_ld:g}_n{nrays}'
    write_adaptive(out, h_ld, nrays, nphi)
    emulate(out, out, training_runs(), max_iter)
    faceon(out)
    return out


def faceon(out):
    """cloudy_rays' face-on map of the emulated run out"""
    from pillardisk.transient.cloudy_rays import faceon_map
    faceon_map(out, 2)


def redo(run, max_iter=400):
    """Re-emulate the dense run <RUNS>/emulated/<run> in place on its own
    plan (current training_runs(): retrains once if they changed) and
    redraw its face-on map"""
    out = os.path.join(RUNS, 'emulated', run)
    emulate(out, out, training_runs(), max_iter)
    faceon(out)


RGB_HS = (0.0, 0.05, 0.2, 0.5, 1.0)     # [ld] rgb: lamp heights
RGB_LINES = (('ha', (1,), r'H\boldmath$\alpha$', 'red'),
             ('mg2', (5, 6), r'Mg\,II', 'lime'),
             ('c4', (3, 4), r'C\,IV', 'dodgerblue'))      # R, G, B


def dense_current(h_ld, nrays=30000, nphi=256, max_iter=400):
    """<RUNS>/emulated/h<h>_n<nrays>, made (dense) if missing and
    re-emulated (redo) if its predictions predate TABLE_VERSION"""
    from pillardisk.transient.cloudy_rays import write_adaptive
    out = f'{RUNS}/emulated/h{h_ld:g}_n{nrays}'
    pred = os.path.join(out, 'emulated_pred.npz')
    ver = -1
    if os.path.exists(pred):
        with np.load(pred) as d:
            ver = int(d['version']) if 'version' in d.files else -1
    if not os.path.exists(os.path.join(out, 'rays_meta.npz')):
        print(f"h = {h_ld:g} ld: no dense run, making {out}", flush=True)
        write_adaptive(out, h_ld, nrays, nphi)
        emulate(out, out, training_runs(), max_iter)
    elif ver != TABLE_VERSION:
        print(f"h = {h_ld:g} ld: {out} predates table version "
              f"{TABLE_VERSION}, re-emulating", flush=True)
        emulate(out, out, training_runs(), max_iter)
    return out


def faceon_cached(out, sim, incl=None, variant=None):
    """cloudy_rays.faceon_surface of RGB_LINES for run out: (S (3, nph,
    nr), ok). incl [deg] given: observed, each ray's line times its escape
    weight toward an observer at incl (line_escape, escape variant
    variant of CIV_VARIANTS). Cached in <out>/faceon_rgb[_i<incl>_<variant>]
    .npz, redone when emulated.npz is newer."""
    from pillardisk.transient.cloudy_rays import (
        faceon_surface, line_escape, CIV_VARIANTS, to_physical)
    tag = '' if incl is None else f'_i{incl:g}_{variant}'
    path = os.path.join(out, f'faceon_rgb{tag}.npz')
    src = os.path.join(out, 'emulated.npz')
    if os.path.exists(path) and os.path.getmtime(path) > \
            os.path.getmtime(src):
        d = np.load(path)
        return d['S'], d['ok']
    rws = [None] * len(RGB_LINES)
    if incl is not None:
        var = [v for v in CIV_VARIANTS if v[0] == variant]
        phys = to_physical(sim, CONFIG)
        rws = []
        for key, _, _, _ in RGB_LINES:
            names, _, ftab = line_escape(out, sim, phys, key, var,
                                         incs=(float(incl),))
            f = ftab[(float(incl), variant)]
            rws.append(dict(zip(names, f)))
            print(f"{os.path.basename(out)} {key}: escape weight "
                  f"median {np.median(f):.2f}", flush=True)
    S, ok, _ = faceon_surface(out, sim, 2, [(c, rw) for (_, c, _, _), rw
                                            in zip(RGB_LINES, rws)])
    np.savez(path, S=S, ok=ok, lines=[k for k, _, _, _ in RGB_LINES])
    return S, ok


def rgb_heights(hs=RGB_HS, nrays=30000, max_iter=400, plot=None,
                incl=None, variant='vtz100'):
    """Face-on RGB composites (R = H alpha, G = Mg II, B = C IV; all the
    emitted light, faceon_surface; or, incl [deg] given, the light
    escaping toward an observer at incl, faceon_cached) of dense emulated
    runs at the lamp
    heights hs [ld], one panel each, one colour stretch for all (each
    channel divided by its largest cell over every panel: hue = line
    balance over RGB_DEX, brightness = strongest line over MAP_DEX).
    Missing or outdated runs are made first (dense_current). Writes
    <RUNS>/emulated/rgb_heights[_i<incl>_<variant>].npz and
    plots/rgb_heights[...].png there."""
    from pillardisk.transient.cloudy_rays import (DUMP, rgb_cells, draw_rgb,
                                                  RGB_DEX, MAP_DEX)
    from pillardisk.transient.sim_spectrum import load_sim
    sim = load_sim(DUMP)
    runs = [dense_current(h, nrays, max_iter=max_iter) for h in hs]
    SS, oks = zip(*[faceon_cached(r, sim, incl, variant) for r in runs])
    tag = '' if incl is None else f'_i{incl:g}_{variant}'
    SS, oks = np.array(SS), np.array(oks)
    x = np.moveaxis(SS, 1, -1)                   # (nh, nph, nr, 3)
    xmax = np.array([np.nanmax(np.where(oks, x[..., k], np.nan))
                     for k in range(3)])
    rgb = rgb_cells(x, xmax)
    plot = plot or os.path.join(RUNS, 'emulated', 'plots',
                                f'rgb_heights{tag}.png')
    os.makedirs(os.path.dirname(plot), exist_ok=True)
    np.savez(os.path.join(os.path.dirname(os.path.dirname(plot)),
                          f'rgb_heights{tag}.npz'), h_ld=np.array(hs), S=SS,
             ok=oks, xmax=xmax, rgb=rgb, runs=np.array(runs),
             lines=[k for k, _, _, _ in RGB_LINES])
    plt = plot_style()
    rf, phf = sim['rf'], sim['phf']
    R, P = np.meshgrid(rf, phf)
    X, Y = R * np.cos(P), R * np.sin(P)
    fig, axs = plt.subplots(1, len(hs), figsize=(6.2 * len(hs), 6.6),
                            sharex=True, sharey=True)
    axs = np.atleast_1d(axs)
    for i, (ax, h) in enumerate(zip(axs, hs)):
        draw_rgb(ax, X, Y, rgb[i], oks[i])
        ax.set_aspect('equal')
        ax.set_xlim(-rf[-1], rf[-1])
        ax.set_ylim(-rf[-1], rf[-1])
        ax.set_title(rf'$\bm{{h = {h:g}}}$\,\textbf{{ld}}', fontsize=24)
        ax.set_xticks([-1, 0, 1])          # none at +-2: panels meet there
        ax.set_xlabel(r'$\bm{x\,[r_0]}$')
        ax.tick_params(which='both', direction='in', labelsize=19)
    axs[0].set_ylabel(r'$\bm{y\,[r_0]}$')
    for j, (_, _, tex, c) in enumerate(RGB_LINES):
        axs[0].text(0.04, 0.96 - 0.08 * j, rf'\textbf{{{tex}}}', color=c,
                    transform=axs[0].transAxes, va='top', fontsize=21)
    axs[-1].text(0.97, 0.96, rf'\textbf{{hue {RGB_DEX:g} dex}}' '\n'
                 rf'\textbf{{brightness {MAP_DEX:g} dex}}', color='0.4',
                 ha='right', va='top', transform=axs[-1].transAxes,
                 fontsize=15)
    if incl is not None:
        vlab = (rf'$\bm{{v_{{t,z}} = {variant[3:]}}}$\,\textbf{{km/s}}'
                if variant.startswith('vtz') else rf'\textbf{{{variant}}}')
        fig.suptitle(rf'\textbf{{observed, $\bm{{i = {incl:g}^\circ}}$, '
                     rf'escape: }}' + vlab, fontsize=24, y=0.98)
    fig.subplots_adjust(wspace=0.04)
    fig.savefig(plot, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {plot}: heights {', '.join(f'{h:g}' for h in hs)} ld")


CMP_LINES = (('ha', (1,), r'H$\bm{\alpha}$'), ('he2', (7,), r'He\,II'),
             ('c4', (3, 4), r'C\,IV'), ('mg2', (5, 6), r'Mg\,II'))
CMP_DEX = 0.5               # [dex] the log-ratio colour range


def compare_faceon(out):
    """Face-on surface brightness of the emulated run out against its
    source run's Cloudy rays (cloudy_rays.faceon_surface on both): rows
    H alpha, He II, C IV, Mg II; columns Cloudy, emulated, log ratio.
    Writes <out>/faceon_compare.npz and <out>/plots/faceon_compare.png"""
    from matplotlib.colors import LogNorm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from pillardisk.transient.cloudy_rays import DUMP, faceon_surface
    from pillardisk.transient.sim_spectrum import load_sim
    src = str(np.load(os.path.join(out, 'emulated.npz'))['src'])
    sim = load_sim(DUMP)
    chans = [(c, None) for _, c, _ in CMP_LINES]
    Sc, okc, _ = faceon_surface(src, sim, 2, chans)
    Se, oke, _ = faceon_surface(out, sim, 2, chans)
    ok = okc & oke
    r0 = CONFIG['R0_LD'] * LD_CM
    rf, phf = sim['rf'], sim['phf']
    area = (sim['r'] * r0 * np.diff(rf) * r0)[None, :] * \
        np.diff(phf)[:, None]
    np.savez(os.path.join(out, 'faceon_compare.npz'), S_cloudy=Sc,
             S_emulated=Se, ok=ok, lines=[k for k, _, _ in CMP_LINES])
    plt = plot_style()
    R, P = np.meshgrid(rf, phf)
    X, Y = R * np.cos(P), R * np.sin(P)
    fig, axs = plt.subplots(len(CMP_LINES), 3, figsize=(21, 6.6 *
                                                         len(CMP_LINES)),
                            sharex=True, sharey=True)
    for i, (key, _, tex) in enumerate(CMP_LINES):
        smax = np.nanmax(np.where(ok, Sc[i], np.nan))
        norm = LogNorm(smax * 1e-5, smax)
        cmap = plt.get_cmap('magma').copy()
        cmap.set_bad('white')
        for j, S in enumerate((Sc[i], Se[i])):
            ax = axs[i, j]
            pc = ax.pcolormesh(X, Y, np.where(ok, np.maximum(S, norm.vmin),
                                              np.nan),
                               norm=norm, cmap=cmap, shading='flat',
                               rasterized=True)
        cax = make_axes_locatable(axs[i, 1]).append_axes('right',
                                                         size='4.5%',
                                                         pad=0.1)
        fig.colorbar(pc, cax=cax).set_label(
            r'$S\,[\rm erg\,s^{-1}\,cm^{-2}]$', fontsize=19)
        make_axes_locatable(axs[i, 0]).append_axes(
            'right', size='4.5%', pad=0.1).set_visible(False)
        lr = np.log10(np.maximum(Se[i], 1e-300) / np.maximum(Sc[i], 1e-300))
        lr = np.where(ok & (Sc[i] > norm.vmin), lr, np.nan)
        pc = axs[i, 2].pcolormesh(X, Y, lr, cmap='RdBu_r', vmin=-CMP_DEX,
                                  vmax=CMP_DEX, shading='flat',
                                  rasterized=True)
        cax = make_axes_locatable(axs[i, 2]).append_axes('right',
                                                         size='4.5%',
                                                         pad=0.1)
        fig.colorbar(pc, cax=cax).set_label(
            r'$\log S_{\rm emu}/S_{\rm Cloudy}$', fontsize=19)
        tot = np.log10(np.nansum(np.where(ok, Se[i], 0) * area)
                       / np.nansum(np.where(ok, Sc[i], 0) * area))
        axs[i, 0].text(0.04, 0.96, r'\textbf{' + tex + r'}', color='k',
                       transform=axs[i, 0].transAxes, va='top',
                       fontsize=26)
        axs[i, 2].text(0.04, 0.96, r'\textbf{total %+.3f dex}' % tot,
                       transform=axs[i, 2].transAxes, va='top',
                       fontsize=19)
        print(f"{key}: emulated / Cloudy total {tot:+.3f} dex; per-cell "
              f"|log ratio| median {np.nanmedian(np.abs(lr)):.3f}, 84% "
              f"{np.nanpercentile(np.abs(lr), 84):.3f} dex")
    for j, t in enumerate((r'\textbf{Cloudy}', r'\textbf{emulated}',
                           r'\textbf{emulated / Cloudy}')):
        axs[0, j].set_title(t, fontsize=26)
    for ax in axs.ravel():
        ax.set_aspect('equal')
        ax.set_xlim(-rf[-1], rf[-1])
        ax.set_ylim(-rf[-1], rf[-1])
        ax.set_yticks([-1, 0, 1, 2])         # none at -2: it meets x = -2
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', labelsize=19,
                       top=True, right=True)
        ax.tick_params(which='major', length=10, width=2)
        ax.tick_params(which='minor', length=5, width=1.5)
    for ax in axs[-1]:
        ax.set_xlabel(r'$x\,[r_0]$')
    for ax in axs[:, 0]:
        ax.set_ylabel(r'$y\,[r_0]$')
    fig.tight_layout(h_pad=0.3, w_pad=0.3)
    pdir = os.path.join(out, 'plots')
    os.makedirs(pdir, exist_ok=True)
    path = os.path.join(pdir, 'faceon_compare.png')
    fig.savefig(path, dpi=100, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('step', choices=['table', 'validate', 'heights',
                                       'heightsplot', 'emulate',
                                       'compare', 'dense', 'synth',
                                       'gapfill', 'redo', 'rgb'])
    ap.add_argument('--out', default=None,
                    help='table / validate: the run directory')
    ap.add_argument('--runs', nargs='+', default=list(HEIGHT_RUNS),
                    help='heights: one run per lamp height')
    ap.add_argument('--src', default=f'{RUNS}/edd10_h0.5_adapt',
                    help='emulate: the run whose ray plan is emulated')
    ap.add_argument('--train', nargs='+', default=None,
                    help='emulate: training runs (default: the heights '
                    'runs but --src)')
    ap.add_argument('--val', choices=sorted(VAL_RUNS), default=None,
                    help='emulate / compare: a Cloudy run to check the '
                    'emulator on (emulated with that height held out, into '
                    '<RUNS>/emulated/<val>_val; emulate then compares)')
    ap.add_argument('--h', type=float, default=0.5,
                    help='dense: lamp height [ld]')
    ap.add_argument('--nrays', type=int, default=30000,
                    help='dense: rays in the plan')
    ap.add_argument('--nphi', type=int, default=256,
                    help='dense: phi columns (at most the sim phi cells)')
    ap.add_argument('--hs', type=float, nargs='+', default=list(RGB_HS),
                    help='rgb: lamp heights [ld]')
    ap.add_argument('--incl', type=float, default=None,
                    help='rgb: observed at this inclination [deg] (escape '
                    'weights; default: all the emitted light)')
    ap.add_argument('--variant', default='vtz100',
                    help='rgb --incl: escape variant (cloudy_rays '
                    'CIV_VARIANTS)')
    ap.add_argument('--run', default=None,
                    help='redo / synth: a run name under <RUNS>/emulated')
    ap.add_argument('--n', type=int, default=GAP_N,
                    help='gapfill: Cloudy rays to add')
    ap.add_argument('--ncpu', type=int, default=GAP_NCPU,
                    help='gapfill: cores for Cloudy')
    ap.add_argument('--norun', action='store_true',
                    help='gapfill: write the decks only')
    ap.add_argument('--smooth', type=int, default=SMOOTH_N,
                    help='synth: running-median rays along a column')
    ap.add_argument('--max-iter', type=int, default=400,
                    help='boosting iterations per target model')
    a = ap.parse_args()
    if a.run is not None and not os.path.isdir(
            os.path.join(RUNS, 'emulated', a.run)):
        have = sorted(os.listdir(os.path.join(RUNS, 'emulated')))
        ap.error(f"no run {a.run!r} under {RUNS}/emulated; have: "
                 + ', '.join(have))
    if a.step == 'table':
        build_table(a.out)
    elif a.step == 'validate':
        validate(a.out, a.max_iter)
    elif a.step == 'emulate':
        src = VAL_RUNS[a.val] if a.val else a.src
        out = a.out or f'{RUNS}/emulated/{a.val or "h0.5"}_val'
        emulate(src, out, a.train, a.max_iter)
        if a.val:
            compare_faceon(out)
    elif a.step == 'redo':
        redo(a.run, a.max_iter)
    elif a.step == 'gapfill':
        gapfill(a.src, a.n, a.out, a.ncpu, not a.norun)
    elif a.step == 'synth':
        out = a.out or os.path.join(RUNS, 'emulated', a.run)
        synth(out, a.smooth)
        faceon(out)
    elif a.step == 'dense':
        dense(a.h, a.nrays, a.nphi, a.out, a.max_iter)
    elif a.step == 'compare':
        compare_faceon(a.out or f'{RUNS}/emulated/{a.val or "h0.5"}_val')
    elif a.step == 'rgb':
        rgb_heights(a.hs, a.nrays, a.max_iter, incl=a.incl,
                    variant=a.variant)
    elif a.step == 'heights':
        validate_heights(a.runs, a.max_iter)
    else:
        heights_plot()


if __name__ == '__main__':
    main()
