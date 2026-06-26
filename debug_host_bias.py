"""Diagnose the host-recovery bias for pillar vs bowl flux-flux.

Imports the already-computed arrays from test_fluxflux_mock. For each disc and
band it prints: the campaign linear fit (a,b), the curvature of the F-X locus
(max deviation of the full locus from the campaign line, in mJy), the anchor
X_non-var, the true lamp-off driver X(s=0), and the recovered host under TWO
anchors -- the F=0 anchor actually used, and the true lamp-off X(s=0).
"""
import numpy as np
import test_fluxflux_mock as m

bands = m.bands


def report(name, F_h, X, fit, anchor, host, floor, true_host):
    xoff = X[0]                                   # driver at s=0 (lamp off)
    xfaint = X[m.camp].min()                      # faintest OBSERVED campaign state
    print(f"\n=== {name} ===  anchor X_nv={anchor:.3f}  lamp-off X(s=0)={xoff:.3f}"
          f"  X_faint={xfaint:.3f}")
    print(f"{'lam':>6} {'a':>7} {'b':>7} {'rec@nv':>7} {'rec@faint':>9} "
          f"{'rec@off':>8} {'true':>7} {'host_in':>7}")
    for k, lam in enumerate(bands):
        a, b = fit[k]
        rec_nv = a * anchor + b
        rec_faint = a * xfaint + b
        rec_off = a * xoff + b
        print(f"{lam:6.0f} {a:7.3f} {b:7.3f} {rec_nv:7.3f} {rec_faint:9.3f} "
              f"{rec_off:8.3f} {true_host[k]:7.3f} {host[k]:7.3f}")


report("PILLAR", m.F_states_h, m.Xh, m.phr, m.Xnv_pr,
       m.HOST, m.F_const_b, m.F_const_h)
report("BOWL", m.F_bowl_h, m.Xhb, m.bhr, m.Xnv_br,
       m.HOST, m.F_const_bowl, m.F_const_bowl + m.HOST)
