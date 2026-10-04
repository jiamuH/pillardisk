"""Check that the light-curve plot runs on the cached single-pillar flux maps
(the step that crashed with np.trapezoid on numpy 1.26). A few seconds.

Run:  python3 debug_lightcurves_from_cache.py
"""
import numpy as np

from pillardisk.pillar_line_time_cloudy import plot_time_evolving_lightcurves

c = np.load('plots/time_evolving_maps_1pillar.npz')
keys = [str(k) for k in c['line_keys']]
plot_time_evolving_lightcurves({k: c['lam_' + k] for k in keys}, c['time'],
                               {k: c['with_' + k] for k in keys},
                               filename='plots/debug_lightcurves_1pillar.png')
