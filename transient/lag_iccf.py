"""Interpolated cross-correlation (ICCF) lag measurement between ZTF bands.

Measures the inter-band continuum lag of the transient quasar with the
standard ICCF method (Peterson et al. 1998): the CCF is computed by
linearly interpolating one light curve onto the shifted time grid of the
other (both directions, averaged), the lag is the centroid of the region
with r > 0.8 r_max, and uncertainties come from flux-randomization /
random-subset-selection (FR/RSS) Monte Carlo.

Sign convention: CCF(tau) = corr[g(t), red(t + tau)], so a positive lag
means the bluer band leads (disk reverberation expectation).

Run: python3 transient/lag_iccf.py
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(HERE, 'data')
PLOTDIR = os.path.join(HERE, 'plots')

TAU_GRID = np.arange(-100.0, 100.0 + 0.25, 0.25)
NMC = 1000
MIN_PAIRS = 20
PAIRS = [('g', 'r'), ('g', 'i')]
WINDOWS = [('flare', 58200.0, 59300.0), ('postflare', 59300.0, np.inf),
           ('full', -np.inf, np.inf)]
# Long-timescale trend removal: Gaussian-weighted running mean with this
# timescale [days] is subtracted in 'detrended' mode, so the CCF is driven
# by the short-timescale variability that carries the reverberation lag.
DETREND_SIGMA = 75.0
MODES = [('raw', None), ('detrended', DETREND_SIGMA)]
RNG = np.random.default_rng(42)


def load_lc(band, mjd_min, mjd_max):
    lc = np.load(os.path.join(DATADIR, f'lc_ztf_{band}.npz'))
    m = (lc['mjd'] >= mjd_min) & (lc['mjd'] <= mjd_max)
    return lc['mjd'][m], lc['flux_mjy'][m], lc['fluxerr_mjy'][m]


def detrend(t, f, sigma):
    """Subtract a Gaussian-weighted running mean with timescale sigma [days]."""
    w = np.exp(-0.5 * ((t[:, None] - t[None, :]) / sigma) ** 2)
    trend = (w @ f) / w.sum(axis=1)
    return f - trend


def ccf_one_direction(t1, f1, t2, f2, tau):
    """Correlate f1(t1) against f2 interpolated at t1 + tau."""
    tshift = t1 + tau
    m = (tshift >= t2[0]) & (tshift <= t2[-1])
    if m.sum() < MIN_PAIRS:
        return np.nan
    a = f1[m]
    b = np.interp(tshift[m], t2, f2)
    sa, sb = a.std(), b.std()
    if sa == 0 or sb == 0:
        return np.nan
    return ((a - a.mean()) * (b - b.mean())).mean() / (sa * sb)


def iccf(t1, f1, t2, f2, tau_grid=TAU_GRID):
    """Two-branch averaged ICCF; branch 2 interpolates LC1 at t2 - tau."""
    r = np.full(len(tau_grid), np.nan)
    for i, tau in enumerate(tau_grid):
        r1 = ccf_one_direction(t1, f1, t2, f2, tau)
        r2 = ccf_one_direction(t2, f2, t1, f1, -tau)
        vals = [v for v in (r1, r2) if np.isfinite(v)]
        if vals:
            r[i] = np.mean(vals)
    return r


def peak_and_centroid(tau_grid, r):
    """Peak lag and the centroid of the contiguous region with r > 0.8 r_max."""
    if not np.any(np.isfinite(r)):
        return np.nan, np.nan, np.nan
    ipk = np.nanargmax(r)
    rmax = r[ipk]
    above = np.isfinite(r) & (r > 0.8 * rmax)
    lo = ipk
    while lo > 0 and above[lo - 1]:
        lo -= 1
    hi = ipk
    while hi < len(r) - 1 and above[hi + 1]:
        hi += 1
    seg_r = r[lo:hi + 1]
    seg_tau = tau_grid[lo:hi + 1]
    centroid = np.sum(seg_tau * seg_r) / np.sum(seg_r)
    return tau_grid[ipk], centroid, rmax


def fr_rss(t, f, e):
    """One FR/RSS realization: subset with replacement + flux randomization."""
    idx = np.unique(RNG.integers(0, len(t), len(t)))
    return t[idx], f[idx] + RNG.normal(0.0, e[idx]), e[idx]


def run_pair(band1, band2, mjd_min, mjd_max, detrend_sigma=None):
    t1, f1, e1 = load_lc(band1, mjd_min, mjd_max)
    t2, f2, e2 = load_lc(band2, mjd_min, mjd_max)
    if detrend_sigma is not None:
        f1 = detrend(t1, f1, detrend_sigma)
        f2 = detrend(t2, f2, detrend_sigma)
    r = iccf(t1, f1, t2, f2)
    tau_peak, tau_cent, rmax = peak_and_centroid(TAU_GRID, r)

    cents, peaks = [], []
    for _ in range(NMC):
        rr = iccf(*fr_rss(t1, f1, e1)[:2], *fr_rss(t2, f2, e2)[:2])
        pk, ct, _ = peak_and_centroid(TAU_GRID, rr)
        if np.isfinite(ct):
            cents.append(ct)
            peaks.append(pk)
    cents, peaks = np.array(cents), np.array(peaks)
    return dict(r=r, tau_peak=tau_peak, tau_cent=tau_cent, rmax=rmax,
                cents=cents, peaks=peaks,
                n1=len(t1), n2=len(t2))


def pct(dist):
    lo, med, hi = np.percentile(dist, [15.87, 50.0, 84.13])
    return med, med - lo, hi - med


def plot_pair(band1, band2, results, mode):
    fig, ax = plt.subplots(figsize=(9, 6.5))
    colors = {'flare': 'orangered', 'postflare': 'seagreen',
              'full': 'royalblue'}
    for wname, res in results.items():
        med, elo, ehi = pct(res['cents'])
        pmed = np.median(res['peaks'])
        label = (rf'$\rm {wname}$: '
                 rf'$\tau_{{\rm cent}} = {med:+.1f}_{{-{elo:.1f}}}'
                 rf'^{{+{ehi:.1f}}}$, '
                 rf'$\tau_{{\rm peak}} = {pmed:+.1f}~\rm d$, '
                 rf'$r_{{\rm max}} = {float(res["rmax"]):.2f}$')
        ax.plot(TAU_GRID, res['r'], lw=3, alpha=0.8, color=colors[wname],
                label=label)
        ax.axvspan(med - elo, med + ehi, color=colors[wname], alpha=0.15, lw=0)
        ax.axvline(med, color=colors[wname], ls='--', lw=1.5)
        hist, edges = np.histogram(res['cents'], bins=40,
                                   range=(TAU_GRID[0], TAU_GRID[-1]))
        if hist.max() > 0:
            centers = 0.5 * (edges[1:] + edges[:-1])
            ax.fill_between(centers, hist / hist.max() * 0.3, 0,
                            color=colors[wname], alpha=0.25, lw=0,
                            step='mid')
    ax.axvline(0.0, color='gray', lw=1, ls=':')
    ax.axhline(0.0, color='gray', lw=1, ls=':')
    ax.set_xlabel(r'$\rm Lag~\tau~[days]$', fontsize=18)
    ax.set_ylabel(rf'$\rm CCF~r(\tau)~~({band1}~vs~{band2},~{mode})$',
                  fontsize=18)
    ax.set_xlim(TAU_GRID[0], TAU_GRID[-1])
    if mode == 'raw':
        ax.set_ylim(bottom=0.0)
    ax.legend(fontsize=12, frameon=False, loc='lower left',
              bbox_to_anchor=(0.0, 1.01))
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    suffix = '' if mode == 'raw' else f'_{mode}'
    out = os.path.join(PLOTDIR, f'transient_iccf_{band1}_{band2}{suffix}.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


def cache_path(band1, band2, mode):
    return os.path.join(DATADIR, f'iccf_cache_{band1}{band2}_{mode}.npz')


CACHE_KEYS = ('r', 'cents', 'peaks', 'rmax', 'n1', 'n2')


def main():
    # --plot-only reuses cached CCFs and Monte-Carlo distributions, so
    # plot tweaks do not redo the (minutes-long) FR/RSS computation.
    plot_only = '--plot-only' in sys.argv
    os.makedirs(PLOTDIR, exist_ok=True)
    print(f'{"pair":>6s} {"mode":>10s} {"window":>7s} {"npts":>10s} '
          f'{"rmax":>6s} {"peak":>16s} {"centroid":>16s}   '
          f'(days; positive = blue leads)')
    for band1, band2 in PAIRS:
        for mode, dsigma in MODES:
            cpath = cache_path(band1, band2, mode)
            if plot_only and os.path.exists(cpath):
                dat = np.load(cpath)
                results = {wname: {k: dat[f'{wname}_{k}']
                                   for k in CACHE_KEYS}
                           for wname, _, _ in WINDOWS}
            else:
                results = {}
                for wname, mjd_min, mjd_max in WINDOWS:
                    results[wname] = run_pair(band1, band2, mjd_min, mjd_max,
                                              detrend_sigma=dsigma)
                np.savez(cpath, **{f'{wname}_{k}': res[k]
                                   for wname, res in results.items()
                                   for k in CACHE_KEYS})
            for wname, res in results.items():
                pmed, plo, phi = pct(res['peaks'])
                cmed, clo, chi = pct(res['cents'])
                print(f'{band1}-{band2:>4s} {mode:>10s} {wname:>7s} '
                      f'{int(res["n1"]):>4d}/{int(res["n2"]):<4d} '
                      f'{float(res["rmax"]):>6.3f} '
                      f'{pmed:>+7.2f} -{plo:.2f}/+{phi:.2f} '
                      f'{cmed:>+7.2f} -{clo:.2f}/+{chi:.2f}')
            plot_pair(band1, band2, results, mode)


if __name__ == '__main__':
    main()
