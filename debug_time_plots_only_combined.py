"""Check that the time-evolving plotting step now makes only the combined
figure (Figures 7 and 8) by default, using the cached single-pillar maps.
A few seconds.

Run:  python3 debug_time_plots_only_combined.py
"""
import glob
import os

import numpy as np

from pillardisk.pillar_line_time_cloudy import _plot_time_evolving_figures

OUT = 'plots/_test_time_plots'
os.makedirs(OUT, exist_ok=True)
config = {'plotting': {k: f'{OUT}/{k}.png' for k in (
    'velocity_time_cloudy_filename', 'velocity_time_cloudy_diff_filename',
    'velocity_time_cloudy_residual_filename', 'velocity_time_cloudy_combined_filename',
    'velocity_time_cloudy_lightcurves_filename')}}
c = np.load('plots/time_evolving_maps_1pillar.npz')
keys = [str(k) for k in c['line_keys']]
_plot_time_evolving_figures(config, {k: c['lam_' + k] for k in keys}, c['time'],
                            {k: c['with_' + k] for k in keys},
                            {k: c['no_' + k] for k in keys},
                            {k: float(c['lambda0_' + k]) for k in keys}, True)
print('files made:', sorted(os.path.basename(f) for f in glob.glob(f'{OUT}/*.png')))
