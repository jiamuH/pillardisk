"""PyROA lag fit to the ZTF g/r/i light curves of the transient quasar.

Prepares PyROA-format data files from transient/data/lc_ztf_{g,r,i}.npz,
runs the joint running-optimal-average MCMC fit (g is the reference band,
so tau_r and tau_i are the lags of r and i behind g, positive = g leads),
prints the lag posteriors, and saves a posterior plot.

PyROA outputs (samples_flat.obj etc.) land in transient/data/pyroa_out/.

Run in the pyroa conda env:
    /Users/jiamuh/miniconda3/envs/pyroa/bin/python3 /Users/jiamuh/python/pillardisk/transient/run_pyroa.py

Options: --nsamples N (default 15000), --nburnin N (default 10000),
         --test (tiny chain, API check only)
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(HERE, 'data')
PYROA_DATADIR = os.path.join(DATADIR, 'pyroa')
OUTDIR = os.path.join(DATADIR, 'pyroa_out')
PLOTDIR = os.path.join(HERE, 'plots')

OBJ = 'transientqso'
FILTERS = ['g', 'r', 'i']
# [A, B multiplicative ranges, tau range [d], delta range [d], extra-var range]
PRIORS = [[0.5, 2.0], [0.5, 2.0], [-50.0, 50.0], [0.01, 10.0], [0.0, 10.0]]


def getarg(name, default):
    if name in sys.argv:
        return int(sys.argv[sys.argv.index(name) + 1])
    return default


def prepare_data():
    os.makedirs(PYROA_DATADIR, exist_ok=True)
    for band in FILTERS:
        lc = np.load(os.path.join(DATADIR, f'lc_ztf_{band}.npz'))
        path = os.path.join(PYROA_DATADIR, f'{OBJ}_{band}.dat')
        np.savetxt(path, np.column_stack(
            [lc['mjd'], lc['flux_mjy'], lc['fluxerr_mjy']]))
        print(f'wrote {path} ({len(lc["mjd"])} points)')


def lag_posteriors(samples_flat):
    """Split the flat chain into per-filter chunks the way PyROA does."""
    ts = np.transpose(samples_flat)
    ts = np.insert(ts, 2, 0.0, axis=0)  # reference-band (g) tau placeholder
    chunks = [ts[i:i + 4] for i in range(0, len(ts), 4)]
    return {band: chunks[i][2] for i, band in enumerate(FILTERS) if i > 0}


def plot_posteriors(taus):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy',
                         'font.size': 20})
    plt.rcParams['text.latex.preamble'] = (r'\usepackage{amsmath} '
                                           r'\usepackage{bm} \boldmath')
    os.makedirs(PLOTDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 6.5))
    colors = {'r': 'orangered', 'i': 'darkgoldenrod'}
    for band, tau in taus.items():
        lo, med, hi = np.percentile(tau, [15.87, 50.0, 84.13])
        label = (rf'$\rm {band}$: $\tau = {med:+.1f}_{{-{med - lo:.1f}}}'
                 rf'^{{+{hi - med:.1f}}}~\rm d$')
        ax.hist(tau, bins=80, density=True, histtype='stepfilled',
                alpha=0.35, color=colors[band])
        ax.hist(tau, bins=80, density=True, histtype='step',
                lw=3, color=colors[band], label=label)
    ax.axvline(0.0, color='gray', lw=1, ls=':')
    ax.set_xlabel(r'$\rm Lag~\tau~[days]~(relative~to~g)$', fontsize=18)
    ax.set_ylabel(r'$\rm Posterior~density$', fontsize=18)
    ax.legend(fontsize=14, frameon=False, loc='upper left')
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    out = os.path.join(PLOTDIR, 'transient_pyroa_lag_posteriors.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


def main():
    test = '--test' in sys.argv
    nsamples = 200 if test else getarg('--nsamples', 15000)
    nburnin = 100 if test else getarg('--nburnin', 10000)

    prepare_data()
    os.makedirs(OUTDIR, exist_ok=True)
    os.chdir(OUTDIR)  # PyROA pickles its outputs into the cwd

    from PyROA import Fit
    # gridsize must be a concrete int: the numba-jitted RunningOptimalAverage
    # cannot type the None default (1000 is PyROA's own fallback value).
    fit = Fit(PYROA_DATADIR + '/', OBJ, FILTERS, PRIORS,
              Nsamples=nsamples, Nburnin=nburnin, add_var=True,
              gridsize=1000)

    taus = lag_posteriors(fit.samples_flat)
    print('\nLag posteriors relative to g (positive = g leads):')
    for band, tau in taus.items():
        lo, med, hi = np.percentile(tau, [15.87, 50.0, 84.13])
        print(f'  {band}: {med:+.2f} -{med - lo:.2f}/+{hi - med:.2f} d')
    np.savez(os.path.join(DATADIR, 'pyroa_lags.npz'),
             **{f'tau_{band}': tau for band, tau in taus.items()})
    plot_posteriors(taus)


if __name__ == '__main__':
    main()
