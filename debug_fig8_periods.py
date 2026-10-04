"""Check the orbital periods quoted for the 100-pillar case (Section 6.3) and
whether the residuals of Figure 8 show a power-spectrum peak near the orbital
period at r_mean = 10 ld. Reads the flux-map cache only; under a second.

Run:  python3 debug_fig8_periods.py
"""
import numpy as np

G, MSUN, C, DAY = 6.674e-8, 1.989e33, 2.998e10, 86400.0
M = 7e7 * MSUN
for r_ld in (4, 5, 10, 20):
    r = r_ld * C * DAY
    T = 2 * np.pi * np.sqrt(r ** 3 / (G * M)) / DAY
    print('T_orb at %2d ld = %6.0f d = %5.1f yr' % (r_ld, T, T / 365.25))

c = np.load('plots/time_evolving_maps_N100.npz')
t = c['time']
print('\nTime grid: %d epochs, %.0f to %.0f d (baseline %.0f d = %.1f yr)'
      % (t.size, t[0], t[-1], t[-1] - t[0], (t[-1] - t[0]) / 365.25))
for k in ('C4', 'Mg2', 'Halpha'):
    f = c['with_' + k]
    res = (f - f.mean(axis=0)) / f.mean()
    p = np.abs(np.fft.rfft(res, axis=0)) ** 2
    freq = np.fft.rfftfreq(t.size, d=t[1] - t[0])
    pw = p.sum(axis=1)[1:]
    per = 1 / freq[1:]
    order = np.argsort(pw)[::-1][:3]
    print('%-6s strongest periods: %s d (fraction of power: %s)'
          % (k, ', '.join('%.0f' % per[i] for i in order),
             ', '.join('%.2f' % (pw[i] / pw.sum()) for i in order)))
