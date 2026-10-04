"""Quick check that the Figure 7/8 combined plot runs with the shared residual
colour range, on small synthetic maps (no Cloudy, a few seconds).

Run:  python3 debug_combined_plot_synthetic.py
"""
import numpy as np

from pillardisk.pillar_line_time_cloudy import plot_time_evolving_raw_residual_combined

t = np.linspace(0, 1200, 40)
l0 = {'C4': 1549.0, 'Mg2': 2798.0, 'Halpha': 6562.8}
lam, flux = {}, {}
for k, amp in (('C4', 1.0), ('Mg2', 0.3), ('Halpha', 0.5)):
    lam[k] = l0[k] * (1 + np.linspace(-0.04, 0.04, 60))
    x = np.linspace(-1, 1, 60)
    stripe = np.exp(-((x[None, :] - 0.6 * np.sin(2 * np.pi * t[:, None] / 800)) / 0.05) ** 2)
    flux[k] = np.exp(-x[None, :] ** 2 / 0.3) + amp * stripe
plot_time_evolving_raw_residual_combined(lam, t, flux, l0,
                                         filename='plots/debug_combined_synthetic.png')
