"""PyROA diagnostic plots for the transient quasar, PG0844-style.

Mirrors ~/python/pg0844/pyroa_pg0844/PG0844_MCMC_plots.py: thin wrappers
around PyROA's built-in plotting utilities (CornerPlot, Lightcurves,
Chains, Convergence) applied to the fit outputs in
transient/data/pyroa_out/. Figures go to transient/plots/.

Usage (pyroa env, from the repo root):
    python3 -m transient.plot_pyroa_diagnostics              # all plots
    python3 -m transient.plot_pyroa_diagnostics corner chains
Plot names: corner, lightcurves, chains, convergence.
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import PyROA

# --- Paths and run parameters ------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(HERE, 'data', 'pyroa') + os.sep
OUTPUTDIR = os.path.join(HERE, 'data', 'pyroa_out') + os.sep
PLOTDIR = os.path.join(HERE, 'plots')

OBJ = 'transientqso'
FILTERS = ['g', 'r', 'i']
DELAY_REF = 'g'
# samples_flat.obj already has the MCMC burn-in discarded by Fit (Nburnin),
# so no further cut is needed here.
BURNIN = 0
BAND_COLORS = ['royalblue', 'orangered', 'darkgoldenrod']

plt.rcParams['savefig.dpi'] = 200


def figname(tag):
    return os.path.join(PLOTDIR, f'transient_pyroa_{tag}.png')


def corner_plots():
    """Corner plot of all parameters, then of the tau (time-delay) subset."""
    PyROA.CornerPlot('all', FILTERS, DELAY_REF, burnin=BURNIN,
                     outputdir=OUTPUTDIR, figname=figname('corner'))
    PyROA.CornerPlot('tau', FILTERS, DELAY_REF, burnin=BURNIN,
                     outputdir=OUTPUTDIR, figname=figname('corner_tau'))


def lightcurves():
    """Light-curve data + best-fit ROA model per band."""
    PyROA.Lightcurves(OBJ, FILTERS, DELAY_REF, datadir=DATADIR,
                      outputdir=OUTPUTDIR, burnin=BURNIN,
                      band_colors=BAND_COLORS, grid=False,
                      show_delay_ref=True, figname=figname('lightcurves'))


def chains():
    """Parameter chain plot from the flattened samples."""
    PyROA.Chains('all', FILTERS, DELAY_REF, burnin=BURNIN,
                 outputdir=OUTPUTDIR, figname=figname('chains'))


def convergence():
    """Autocorrelation-based convergence check."""
    # PyROA hardcodes savefig('pyroa_convergence.pdf') into the cwd.
    cwd = os.getcwd()
    os.chdir(PLOTDIR)
    try:
        PyROA.Convergence(outputdir=OUTPUTDIR)
    finally:
        os.chdir(cwd)
    os.replace(os.path.join(PLOTDIR, 'pyroa_convergence.pdf'),
               os.path.join(PLOTDIR, 'transient_pyroa_convergence.pdf'))


# --- Command-line entry point -------------------------------------------------
PLOTS = {
    'corner': corner_plots,
    'lightcurves': lightcurves,
    'chains': chains,
    'convergence': convergence,
}


def main():
    parser = argparse.ArgumentParser(
        description='Make PyROA post-processing plots for the transient quasar.')
    parser.add_argument(
        'plots', nargs='*', default=list(PLOTS),
        help="plots to make: one or more of {%s}. Default: all."
             % ', '.join(PLOTS))
    args = parser.parse_args()

    os.makedirs(PLOTDIR, exist_ok=True)
    for name in args.plots:
        if name not in PLOTS:
            parser.error("unknown plot '%s'; choose from %s"
                         % (name, ', '.join(PLOTS)))
        print(f'--- {name}')
        PLOTS[name]()


if __name__ == '__main__':
    main()
